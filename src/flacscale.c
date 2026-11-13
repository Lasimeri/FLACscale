/*
 * flacscale.c -- the command: one FLAC, WAV or MP3 in, one upscaled 32-bit FLAC out.
 *
 * Purpose : decode -> repair clipped peaks (Opus LPC model) -> resample to the
 *           target rate in double precision -> 32-bit samples -> FLAC, tags
 *           and cover art carried over.
 *           A lossy input (MP3) is also written out as a plain lossless FLAC
 *           (24 bit, source rate) beside the output; the upscaler itself is fed
 *           from the decoder's floating-point output, not from that file.
 * Passes  : 0  scan    decode once, decide per channel whether it is clipped
 *                      and at what level (skipped for lossy input and with -D)
 *           1  measure run the signal path, find the true peak at the output
 *                      rate (skipped with a fixed -g DB)
 *           2  write   run it again with the smallest uniform attenuation that
 *                      avoids clipping (none when the peak fits) and encode
 * Usage   : flacscale IN [-o OUT.flac] [-r RATE] [-g auto|DB] [-t THREADS] [-n] [-D] [-X]
 */
#include "fu_cliplevel.h"
#include "fu_cutoff.h"
#include "fu_engine.h"
#include "fu_flac.h"
#include "fu_input.h"
#include "fu_pcm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

enum { OUT_BITS = 32, COPY_BITS = 24 };

typedef struct options {
    const char *in, *out;
    uint32_t    rate, threads;
    int         gain_auto, no_copy, no_declip, no_extend;
    double      gain_db;
} options;

typedef struct pass_stats {
    uint64_t frames_in, frames_out, clipped_out, clipped_src;
    double   peak_out, peak_src;
    double   s_decode, s_dsp, s_encode;
    fu_declip_stats declip;
    fu_extend_stats extend;
    int             extending;
} pass_stats;

/* Where one pass sends its samples; a NULL encoder means "measure only". */
typedef struct sinks {
    fu_flac_out *out, *copy;
    double       gain_out, gain_copy;
} sinks;

typedef struct work {
    double  *sig[FU_MAX_CH];
    int32_t *pcm[FU_MAX_CH];
    size_t   cap;
} work;

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static double db(double ratio)
{
    return 20.0 * log10(ratio);
}

static void usage(FILE *f)
{
    fputs("usage: flacscale IN [-o OUT.flac] [-r RATE] [-g auto|DB] [-t THREADS] [-n] [-D] [-X]\n"
          "  IN            a FLAC, WAV, AIFF or MP3 file (recognised by content)\n"
          "  -o OUT.flac   output file            (default: IN's name + .RATEk32.flac, here)\n"
          "  -r RATE       output sample rate, Hz (default: 192000)\n"
          "  -g auto|DB    output gain            (default: auto = only as much attenuation\n"
          "                                        as the true peak needs; DB = fixed, e.g. -3)\n"
          "  -t THREADS    FLAC encoder threads   (default: 4)\n"
          "  -n            lossy input: do not also write the plain 24-bit FLAC conversion\n"
          "  -D            do not repair clipped peaks\n"
          "  -X            do not extend the band above the source's real bandwidth\n",
          f);
}

static int parse(int argc, char **argv, options *o)
{
    *o = (options){ .rate = 192000, .threads = 4, .gain_auto = 1 };
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            usage(stdout);
            exit(0);
        }
        if (!strcmp(a, "-n")) {
            o->no_copy = 1;
        } else if (!strcmp(a, "-X")) {
            o->no_extend = 1;
        } else if (!strcmp(a, "-D")) {
            o->no_declip = 1;
        } else if (a[0] == '-' && a[1] && strchr("orgt", a[1]) && !a[2]) {
            if (++i >= argc) {
                fprintf(stderr, "flacscale: %s needs a value\n", a);
                return -1;
            }
            const char *v = argv[i];
            switch (a[1]) {
            case 'o': o->out = v; break;
            case 'r': o->rate = (uint32_t)strtoul(v, NULL, 10); break;
            case 't': o->threads = (uint32_t)strtoul(v, NULL, 10); break;
            default:
                o->gain_auto = !strcmp(v, "auto");
                o->gain_db = o->gain_auto ? 0.0 : strtod(v, NULL);
                break;
            }
        } else if (a[0] == '-' && a[1]) {
            fprintf(stderr, "flacscale: unknown option %s\n", a);
            return -1;
        } else if (!o->in) {
            o->in = a;
        } else {
            fprintf(stderr, "flacscale: one input file at a time (extra argument: %s)\n", a);
            return -1;
        }
    }
    if (!o->in) {
        fprintf(stderr, "flacscale: no input file\n");
        return -1;
    }
    return 0;
}

/* Length of IN's file name without directory and without its extension
 * (a final ".xxx" of at most 5 characters: .flac, .mp3, .wav, .aiff). */
static size_t stem_of(const char *path, const char **base)
{
    const char *b = strrchr(path, '/');
    b = b ? b + 1 : path;
    const size_t n = strlen(b);
    const char *dot = strrchr(b, '.');
    *base = b;
    return dot && dot != b && strlen(dot) <= 6 ? (size_t)(dot - b) : n;
}

static int work_reserve(work *w, uint32_t nch, size_t frames)
{
    if (frames <= w->cap)
        return 0;
    for (uint32_t ch = 0; ch < nch; ch++) {
        double *s = realloc(w->sig[ch], frames * sizeof *s);
        if (!s)
            return -1;
        w->sig[ch] = s;
        int32_t *p = realloc(w->pcm[ch], frames * sizeof *p);
        if (!p)
            return -1;
        w->pcm[ch] = p;
    }
    w->cap = frames;
    return 0;
}

static void work_free(work *w)
{
    for (int ch = 0; ch < FU_MAX_CH; ch++) {
        free(w->sig[ch]);
        free(w->pcm[ch]);
    }
}

/* Pass 0: decode once; decide each channel's ceilings and where the source's
 * real bandwidth ends. *wall = 0 when any channel reaches Nyquist without a
 * wall, else the highest wall among the channels; *drop is its depth in dB. */
static int scan_source(const options *o, uint32_t nch, uint32_t rate, fu_clipverdict *verdict,
                       double *wall, double *drop)
{
    char err[512];
    fu_format src;
    fu_cliplevel *lv = malloc(nch * sizeof *lv);
    fu_cutoff *cut[FU_MAX_CH] = { 0 };
    fu_input *in = fu_input_open(o->in, &src, err, sizeof err);
    int ready = lv && in;
    for (uint32_t ch = 0; ch < nch; ch++)
        ready &= (cut[ch] = fu_cutoff_new(rate)) != NULL;
    long n = -1;
    if (ready) {
        for (uint32_t ch = 0; ch < nch; ch++)
            fu_cliplevel_reset(&lv[ch]);
        const double *planes[FU_MAX_CH];
        while ((n = fu_input_read(in, planes)) > 0)
            for (uint32_t ch = 0; ch < nch; ch++) {
                fu_cliplevel_feed(&lv[ch], planes[ch], (size_t)n);
                fu_cutoff_feed(cut[ch], planes[ch], (size_t)n);
            }
        int open_band = 0;
        *wall = 0.0;
        for (uint32_t ch = 0; ch < nch; ch++) {
            double d = 0.0;
            const double w = fu_cutoff_decide(cut[ch], &d);
            verdict[ch] = fu_cliplevel_decide(&lv[ch]);
            open_band |= w == 0.0;
            if (w > *wall) {
                *wall = w;
                *drop = d;
            }
        }
        if (open_band)
            *wall = 0.0;
    }
    for (uint32_t ch = 0; ch < nch; ch++)
        fu_cutoff_free(cut[ch]);
    fu_input_close(in);
    free(lv);
    if (n != 0)
        fprintf(stderr, "flacscale: %s: cannot scan the file\n", o->in);
    return n == 0 ? 0 : -1;
}

/* Measure one block; when an encoder is attached, convert and write it. */
static int emit(const double *const *sig, int32_t *const *pcm, uint32_t nch, size_t n,
                fu_flac_out *enc, double gain, unsigned bits, double *peak, uint64_t *clipped,
                pass_stats *st)
{
    double t = now_s();
    for (uint32_t ch = 0; ch < nch; ch++) {
        *peak = fmax(*peak, fu_pcm_peak(sig[ch], n));
        if (enc)
            *clipped += fu_pcm_to_int(sig[ch], pcm[ch], n, gain, bits);
    }
    st->s_dsp += now_s() - t;
    if (!enc)
        return 0;
    t = now_s();
    const int rc = fu_flac_write(enc, (const int32_t *const *)pcm, n);
    st->s_encode += now_s() - t;
    return rc;
}

/* Passes 1 and 2: one trip through the signal path. */
static int run_pass(const options *o, const fu_engine_cfg *cfg, const sinks *k, pass_stats *st)
{
    char err[512];
    fu_format src;
    work w = { 0 };
    int rc = -1;

    memset(st, 0, sizeof *st);
    fu_input *in = fu_input_open(o->in, &src, err, sizeof err);
    if (!in) {
        fprintf(stderr, "flacscale: %s\n", err);
        return -1;
    }
    fu_engine *eng = fu_engine_new(cfg);
    if (!eng) {
        fprintf(stderr, "flacscale: cannot build a %u -> %u Hz signal path; -r must be at least the "
                        "input rate (%u)\n", src.rate, o->rate, src.rate);
        fu_input_close(in);
        return -1;
    }
    for (;;) {
        const double *planes[FU_MAX_CH];
        double t = now_s();
        const long n = fu_input_read(in, planes);
        st->s_decode += now_s() - t;
        if (n < 0) {
            fprintf(stderr, "flacscale: %s: decode error after %llu frames; the file is damaged\n",
                    o->in, (unsigned long long)st->frames_in);
            goto done;
        }
        const size_t need = fu_engine_bound(eng, (size_t)n);
        if (work_reserve(&w, src.channels, need > (size_t)n ? need : (size_t)n)) {
            fprintf(stderr, "flacscale: out of memory\n");
            goto done;
        }
        /* The plain conversion leaves first: it reuses the integer buffers. */
        if (n && emit(planes, w.pcm, src.channels, (size_t)n, k->copy, k->gain_copy, COPY_BITS,
                      &st->peak_src, &st->clipped_src, st))
            goto write_failed;
        /* Self-check: a sane decode never exceeds full scale by 24 dB. NaN fails
         * the comparison too, so garbage cannot be "normalised" into the output. */
        if (!(st->peak_src < 16.0)) {
            fprintf(stderr, "flacscale: %s: the decoder produced invalid samples (peak %g) near "
                            "frame %llu; refusing to write them\n", o->in, st->peak_src,
                    (unsigned long long)st->frames_in);
            goto done;
        }
        t = now_s();
        const size_t m = n ? fu_engine_run(eng, planes, (size_t)n, w.sig, w.cap)
                           : fu_engine_flush(eng, w.sig, w.cap);
        st->s_dsp += now_s() - t;
        st->frames_in += (uint64_t)n;
        st->frames_out += m;
        if (emit((const double *const *)w.sig, w.pcm, src.channels, m, k->out, k->gain_out,
                 OUT_BITS, &st->peak_out, &st->clipped_out, st))
            goto write_failed;
        if (n == 0)
            break;
    }
    st->declip = *fu_engine_declip_stats(eng);
    st->extend = *fu_engine_extend_stats(eng);
    st->extending = fu_engine_extending(eng);
    rc = 0;
    goto done;
write_failed:
    fprintf(stderr, "flacscale: write failed (disk full?)\n");
done:
    fu_engine_free(eng);
    fu_input_close(in);
    work_free(&w);
    return rc;
}

static double file_mib(const char *path)
{
    struct stat sb;
    return stat(path, &sb) == 0 ? (double)sb.st_size / 1048576.0 : 0.0;
}

int main(int argc, char **argv)
{
    options o;
    char outbuf[4096], copybuf[4096], err[512], tag[256];

    if (parse(argc, argv, &o)) {
        usage(stderr);
        return 2;
    }
    const char *base;
    const size_t stem = stem_of(o.in, &base);
    if (!o.out) {
        snprintf(outbuf, sizeof outbuf, "%.*s.%uk%d.flac", (int)stem, base, o.rate / 1000, OUT_BITS);
        o.out = outbuf;
    }

    fu_format src;
    fu_input *probe = fu_input_open(o.in, &src, err, sizeof err);
    if (!probe) {
        fprintf(stderr, "flacscale: %s\n", err);
        return 1;
    }
    const char *kind = fu_input_kind(probe);
    const int lossy = (int)src.lossy, want_copy = lossy && !o.no_copy;
    const int declip = !lossy && !o.no_declip;
    char depth[16];
    snprintf(depth, sizeof depth, lossy ? "lossy" : src.bps ? "%u bit" : "float", src.bps);
    if (want_copy) {
        const char *slash = strrchr(o.out, '/');
        snprintf(copybuf, sizeof copybuf, "%.*s%.*s.flac", slash ? (int)(slash - o.out + 1) : 0,
                 o.out, (int)stem, base);
    }
    const double t_start = now_s();

    /* Pass 0. */
    fu_engine_cfg cfg = fu_engine_cfg_plain(src.rate, o.rate, src.channels);
    fu_clipverdict verdict[FU_MAX_CH] = { 0 };
    double wall = 0.0, wall_drop = 0.0;
    if (scan_source(&o, src.channels, src.rate, verdict, &wall, &wall_drop)) {
        fu_input_close(probe);
        return 1;
    }
    for (uint32_t ch = 0; declip && ch < src.channels; ch++) {
        cfg.hi[ch] = verdict[ch].hi;
        cfg.lo[ch] = verdict[ch].lo;
    }
    if (!o.no_extend) {
        /* With a wall, real content ends there. Without one it runs to the source
         * Nyquist, but only the resampler's flat passband (0.907 of it) is trusted. */
        cfg.src_top = wall > 0.0 ? wall : 0.5 * src.rate;
        cfg.ref_top = (wall > 0.0 ? wall : 0.907 * 0.5 * src.rate) - 250.0;
    }

    /* Pass 1. */
    pass_stats a = { 0 }, b;
    sinks k = { .gain_out = pow(10.0, o.gain_db / 20.0), .gain_copy = 1.0 };
    if (o.gain_auto) {
        const sinks measure = { .gain_out = 1.0, .gain_copy = 1.0 };
        if (run_pass(&o, &cfg, &measure, &a)) {
            fu_input_close(probe);
            return 1;
        }
        const double lim_out = fu_pcm_limit(OUT_BITS), lim_copy = fu_pcm_limit(COPY_BITS);
        k.gain_out  = a.peak_out > lim_out ? lim_out / a.peak_out : 1.0;
        k.gain_copy = a.peak_src > lim_copy ? lim_copy / a.peak_src : 1.0;
    }
    const double t_analyse = now_s() - t_start;

    /* Pass 2. */
    const fu_format dst = {
        .rate = o.rate, .channels = src.channels, .bps = OUT_BITS,
        .frames = (src.frames * o.rate + src.rate - 1) / src.rate,
    };
    const char *const extra[] = { tag };
    snprintf(tag, sizeof tag, "FLACSCALE=%s %u Hz %s -> %u Hz %d bit, gain %+.4f dB", kind, src.rate,
             depth, dst.rate, OUT_BITS, db(k.gain_out));
    k.out = fu_flac_create(o.out, &dst, fu_input_tags(probe), extra, 1, o.threads, err, sizeof err);
    if (k.out && want_copy) {
        const fu_format cfmt = { .rate = src.rate, .channels = src.channels, .bps = COPY_BITS,
                                 .frames = src.frames };
        snprintf(tag, sizeof tag, "FLACSCALE=decoded from %s to %d bit, gain %+.4f dB", kind,
                 COPY_BITS, db(k.gain_copy));
        k.copy = fu_flac_create(copybuf, &cfmt, fu_input_tags(probe), extra, 1, o.threads, err,
                                sizeof err);
    }
    fu_input_close(probe);
    if (!k.out || (want_copy && !k.copy)) {
        fprintf(stderr, "flacscale: %s\n", err);
        return 1;
    }
    int rc = run_pass(&o, &cfg, &k, &b);
    rc |= fu_flac_finish(k.out);
    if (k.copy)
        rc |= fu_flac_finish(k.copy);
    if (rc) {
        fprintf(stderr, "flacscale: %s: output is incomplete\n", o.out);
        return 1;
    }
    const double t_total = now_s() - t_start, dur = (double)b.frames_in / src.rate;

    printf("flacscale: %s\n", o.in);
    printf("  input      %-4s %6u Hz  %-6s  %u ch  %llu frames  %.3f s  %.1f MiB\n", kind, src.rate,
           depth, src.channels, (unsigned long long)b.frames_in, dur, file_mib(o.in));
    if (k.copy)
        printf("  converted  FLAC %6u Hz  %d bit  %u ch  %llu frames  %.1f MiB  peak %+.3f dBFS  "
               "gain %+.4f dB  clipped %llu  -> %s\n", src.rate, COPY_BITS, src.channels,
               (unsigned long long)b.frames_in, file_mib(copybuf), db(b.peak_src), db(k.gain_copy),
               (unsigned long long)b.clipped_src, copybuf);
    for (uint32_t ch = 0; declip && ch < src.channels; ch++) {
        const fu_clipverdict *v = &verdict[ch];
        printf("  clipping   ch %u: %llu samples at the top level (%.0f on the levels inside it), "
               "%llu at the bottom (%.0f): %s\n", ch, (unsigned long long)v->at_hi, v->inside_hi,
               (unsigned long long)v->at_lo, v->inside_lo,
               isfinite(v->hi) || isfinite(v->lo) ? "CLIPPED, repairing" : "clean");
    }
    if (declip && b.declip.pinned)
        printf("  repair     %llu pinned samples in %llu clusters: %llu raised past the ceiling, "
               "%llu held at it, %llu clusters refused; highest repaired peak %+.3f dBFS\n",
               (unsigned long long)b.declip.pinned, (unsigned long long)b.declip.clusters,
               (unsigned long long)b.declip.raised, (unsigned long long)b.declip.held,
               (unsigned long long)b.declip.refused, db(b.declip.peak));
    if (wall > 0.0)
        printf("  bandwidth  the source's content ends at %.0f Hz (a %.1f dB wall)\n", wall, wall_drop);
    else
        printf("  bandwidth  the source's content reaches its Nyquist (%.0f Hz), no wall\n",
               0.5 * src.rate);
    if (b.extending)
        printf("  extension  folded %.0f .. 25000 Hz from the band below %.0f Hz; above 25000 Hz "
               "harmonics only: %llu tonal peaks -> %llu partials (%.1f per second)\n",
               cfg.src_top, cfg.ref_top, (unsigned long long)b.extend.peaks,
               (unsigned long long)b.extend.lobes,
               (double)b.extend.lobes / dur);
    printf("  output     FLAC %6u Hz  %d bit  %u ch  %llu frames  %.1f MiB  -> %s\n", dst.rate,
           OUT_BITS, dst.channels, (unsigned long long)b.frames_out, file_mib(o.out), o.out);
    printf("  true peak  %+.3f dBFS at %u Hz   gain %+.4f dB (%s)   clipped samples %llu\n",
           db(b.peak_out), dst.rate, db(k.gain_out), o.gain_auto ? "auto" : "fixed",
           (unsigned long long)b.clipped_out);
    printf("  time       scan+measure %.2f s | decode %.2f s  dsp %.2f s  encode %.2f s | total "
           "%.2f s  (%.1fx real time)\n", t_analyse, b.s_decode, b.s_dsp, b.s_encode, t_total,
           dur / t_total);
    return 0;
}
