/*
 * dequantcheck.c -- gate: score the dequantizer against a known truth.
 *
 * Purpose : take audio that is finer than the target bit depth (a real 24-bit
 *           master, or a float decode), round an excerpt of it to BITS, run
 *           fu_dequant on the rounded copy, and measure how close each version
 *           is to the original:
 *             SDR = 10*log10( sum(truth^2) / sum((truth - estimate)^2) )
 *           "in" scores the rounded signal, "out" the dequantized one.
 *           -a attenuates the truth first: a passage at -40 dBFS in a 16-bit
 *           file has about 9 bits to live on, and that is where rounding hurts.
 *           -d adds TPDF dither before rounding (as mastering usually does) and
 *           widens the consistency interval to the 1.5 steps that implies.
 * Input   : any file fu_input reads. If the truth already sits on the BITS
 *           grid, there is nothing to recover and the tool says so.
 * Usage   : dequantcheck FILE [-b BITS] [-a ATTEN_DB] [-d] [-s START_S] [-l LEN_S]
 *                             [-n PASSES] [-w HALF] [-t THREADS]
 *           defaults: 16 bits, 0 dB, 20 s from the middle, library settings
 */
#include "fu_dequant.h"
#include "fu_input.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct options {
    const char *file;
    int         bits, dither;
    double      atten_db, start_s, len_s, half;
    fu_dequant_cfg cfg;
} options;

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int parse(int argc, char **argv, options *o)
{
    *o = (options){ .bits = 16, .start_s = -1.0, .len_s = 20.0, .half = -1.0,
                    .cfg = fu_dequant_cfg_default() };
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (a[0] != '-' || !a[1]) {
            if (o->file)
                return -1;
            o->file = a;
        } else if (a[1] == 'd' && !a[2]) {
            o->dither = 1;
        } else if (strchr("baslnwt", a[1]) && !a[2] && i + 1 < argc) {
            const char *v = argv[++i];
            switch (a[1]) {
            case 'b': o->bits = atoi(v); break;
            case 'a': o->atten_db = strtod(v, NULL); break;
            case 's': o->start_s = strtod(v, NULL); break;
            case 'l': o->len_s = strtod(v, NULL); break;
            case 'n': o->cfg.passes = atoi(v); break;
            case 'w': o->half = strtod(v, NULL); break;
            default:  o->cfg.threads = atoi(v); break;
            }
        } else {
            return -1;
        }
    }
    return o->file && o->bits >= 2 && o->bits <= 24 && o->len_s > 0.0 ? 0 : -1;
}

/* Read frames [first, first + count) of every channel; returns the frames read. */
static size_t read_excerpt(fu_input *in, uint32_t nch, size_t first, size_t count, double **ch)
{
    size_t pos = 0, got_total = 0;
    for (;;) {
        const double *planes[FU_MAX_CH];
        const long got = fu_input_read(in, planes);
        if (got <= 0)
            break;
        for (long i = 0; i < got; i++, pos++) {
            if (pos < first || got_total >= count)
                continue;
            for (uint32_t c = 0; c < nch; c++)
                ch[c][got_total] = planes[c][i];
            got_total++;
        }
        if (got_total >= count)
            break;
    }
    return got_total;
}

int main(int argc, char **argv)
{
    options o;
    if (parse(argc, argv, &o)) {
        fprintf(stderr, "usage: dequantcheck FILE [-b BITS] [-a ATTEN_DB] [-d] [-s START_S] "
                        "[-l LEN_S] [-n PASSES] [-w HALF] [-t THREADS]\n");
        return 2;
    }
    char err[512];
    fu_format fmt;
    fu_input *in = fu_input_open(o.file, &fmt, err, sizeof err);
    if (!in) {
        fprintf(stderr, "dequantcheck: %s\n", err);
        return 1;
    }
    const double total_s = (double)fmt.frames / fmt.rate;
    if (o.len_s > total_s)
        o.len_s = total_s;
    if (o.start_s < 0.0)
        o.start_s = 0.5 * (total_s - o.len_s);
    const size_t first = (size_t)(o.start_s * fmt.rate), want = (size_t)(o.len_s * fmt.rate);
    double *truth[FU_MAX_CH] = { 0 }, *est[FU_MAX_CH] = { 0 };
    for (uint32_t c = 0; c < fmt.channels; c++) {
        truth[c] = malloc(want * sizeof(double));
        est[c] = malloc(want * sizeof(double));
        if (!truth[c] || !est[c]) {
            fprintf(stderr, "dequantcheck: out of memory\n");
            return 1;
        }
    }
    const size_t n = read_excerpt(in, fmt.channels, first, want, truth);
    fu_input_close(in);
    if (n == 0) {
        fprintf(stderr, "dequantcheck: %s: no frames in the excerpt\n", o.file);
        return 1;
    }

    const double step = ldexp(1.0, 1 - o.bits), gain = pow(10.0, o.atten_db / 20.0);
    const double half = o.half > 0.0 ? o.half : o.dither ? 1.5 : 0.5;
    o.cfg.half = half;
    o.cfg.noise = o.dither ? 0.25 : 1.0 / 12.0;
    uint64_t r = 0x9E3779B97F4A7C15ull;
    double tt = 0.0, e_in = 0.0, e_out = 0.0, worst = 0.0, secs = 0.0;
    size_t on_grid = 0;
    for (uint32_t c = 0; c < fmt.channels; c++) {
        double *t = truth[c], *q = est[c];
        for (size_t i = 0; i < n; i++) {
            t[i] *= gain;
            double d = 0.0;
            if (o.dither)
                for (int k = 0; k < 2; k++) {       /* TPDF: two uniform values, +-0.5 step each */
                    r ^= r << 13;
                    r ^= r >> 7;
                    r ^= r << 17;
                    d += (double)(r >> 11) * 0x1p-53 - 0.5;
                }
            const double level = fmin(fmax(nearbyint(t[i] / step + d), -1.0 / step), 1.0 / step - 1.0);
            q[i] = level * step;
            on_grid += t[i] == q[i];
            tt += t[i] * t[i];
            e_in += (t[i] - q[i]) * (t[i] - q[i]);
        }
        double *work = malloc(n * sizeof *work);
        if (!work)
            return 1;
        memcpy(work, q, n * sizeof *work);
        const double t0 = now_s();
        if (fu_dequant(work, n, step, &o.cfg)) {
            fprintf(stderr, "dequantcheck: out of memory in the dequantizer\n");
            return 1;
        }
        secs += now_s() - t0;
        for (size_t i = 0; i < n; i++) {
            e_out += (t[i] - work[i]) * (t[i] - work[i]);
            worst = fmax(worst, fabs(work[i] - q[i]) / step);
        }
        free(work);
    }

    const double sdr_in = 10.0 * log10(tt / (e_in + 1e-300));
    const double sdr_out = 10.0 * log10(tt / (e_out + 1e-300));
    const double rms_db = 10.0 * log10(tt / ((double)n * fmt.channels) + 1e-300);
    printf("dequantcheck: %u Hz, %u ch, %.3f s from %.3f s; level %+.1f dB (RMS %.1f dBFS); "
           "%d bit%s; box +-%.2f steps; %d pass(es)\n",
           fmt.rate, fmt.channels, (double)n / fmt.rate, o.start_s, o.atten_db, rms_db, o.bits,
           o.dither ? " + TPDF dither" : "", half, o.cfg.passes);
    if (on_grid == n * fmt.channels) {
        printf("  the truth already sits on the %d-bit grid: nothing to recover, no score\n", o.bits);
        return 0;
    }
    printf("  SDR rounded %.2f dB   SDR dequantized %.2f dB   gain %+.2f dB   "
           "largest move %.3f steps   %.2f s (%.1fx real time)\n",
           sdr_in, sdr_out, sdr_out - sdr_in, worst, secs,
           (double)n * fmt.channels / fmt.rate / secs);
    return 0;
}
