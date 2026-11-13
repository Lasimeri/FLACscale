/*
 * fu_extend.c -- band extension (see fu_extend.h).
 *
 * Transform: frames of N = 4096 samples, hop H = N/4, periodic Hann on analysis
 *            and synthesis (sum of squares over the 4 overlaps = 1.5).
 * Phases   : bins are multiplied by (-1)^k, which moves the phase reference to
 *            the centre of the frame. A stationary partial then has one phase
 *            across its whole main lobe, a bin copy by D bins is coherent
 *            between frames when D is a multiple of N/H = 4, and a harmonic
 *            h gets phase h * phi with no window term to correct.
 * Lobe     : centred Hann kernel, K(d) = sinc(d) / (1 - d^2); 16 bins are written.
 * Guard    : per frame, synthesised power above fold_top <= -30 dB of the real
 *            signal power (a wrong guess must not be able to load a tweeter).
 * Stream   : input history N - H frames; after each hop one frame is analysed,
 *            its synthesis is overlap-added, and the oldest H frames leave as
 *            dry + extension. The first N - H frames out are the start-up
 *            delay and are dropped; flush feeds zeros to push the tail out.
 */
#include "fu_extend.h"
#include "fu_audio.h"

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { N = 4096, H = N / 4, NBIN = N / 2 + 1, SMOOTH = 8, FLOOR_HALF = 24, FLOOR_SKIP = 3, LOBE = 8 };
#define TONAL_RATIO 31.6        /* 15 dB above the local floor, as a power ratio */
#define OLA_NORM    (1.0 / 1.5)
#define ULTRA_CAP   1e-3        /* -30 dB: ceiling on synthesised power above fold_top */

typedef struct chan {
    double *in;                 /* N - H history + new frames */
    double *ola;                /* N frames of overlap-add */
    size_t  fill;               /* frames in `in` */
    double  slope;              /* smoothed envelope slope, dB/kHz */
    double  dens_lo, dens_hi;   /* running densities of the two octaves below src_top */
    int     primed;
} chan;

struct fu_extend {
    fu_extend_cfg   cfg;
    fu_extend_stats st;
    chan            ch[FU_MAX_CH];
    size_t          cap;                /* frames allocated in ch[].in */
    uint64_t        total_in, total_out;
    size_t          skip;               /* start-up frames still to drop */
    double          hz;                 /* Hz per bin */
    int             k_top;              /* first bin above the usable source band */
    int             k_src, k_fold, k_fold_end, src_bins;    /* fold geometry */
    int             k_harm, k_max, k_peak0;                 /* harmonic geometry */
    double          fold_gain, harm_gain;
    double          win[N], pw[NBIN], cs[NBIN + 1];
    double         *frame;
    fftw_complex   *spec, *synth;
    fftw_plan       fwd, inv;
};

fu_extend_cfg fu_extend_cfg_default(uint32_t rate, uint32_t channels, double src_top,
                                    double ref_top)
{
    return (fu_extend_cfg){ .rate = rate, .channels = channels, .src_top = src_top, .ref_top = ref_top,
                            .fold_top = 25000.0, .max_harm = 6 };
}

static double kernel(double d)
{
    const double a = fabs(d);
    if (a < 1e-9)
        return 1.0;
    if (fabs(a - 1.0) < 1e-9)
        return 0.5;
    return sin(M_PI * d) / (M_PI * d * (1.0 - d * d));
}

fu_extend *fu_extend_new(const fu_extend_cfg *cfg)
{
    fu_extend *x = calloc(1, sizeof *x);
    if (!x)
        return NULL;
    x->cfg = *cfg;
    x->hz = (double)cfg->rate / N;
    x->k_top = (int)floor(cfg->ref_top / x->hz);
    x->k_fold_end = (int)floor(cfg->fold_top / x->hz);
    x->k_max = (int)floor(0.47 * cfg->rate / x->hz);
    x->k_harm = x->k_fold_end;
    x->k_peak0 = x->k_harm / cfg->max_harm;

    /* Fold geometry: a source band of src_bins (multiple of 4) ending at k_top
     * is tiled upwards from k_fold; k_fold - k_src is a multiple of 4. */
    const int k_floor = (int)ceil(4000.0 / x->hz);
    int k_fold = (int)ceil(cfg->src_top / x->hz);
    int want = x->k_fold_end - k_fold;
    int have = x->k_top - k_floor;
    x->src_bins = (want < have ? want : have) & ~3;
    x->k_src = x->k_top - x->src_bins;
    k_fold += (4 - ((k_fold - x->k_src) & 3)) & 3;
    x->k_fold = k_fold;

    if (cfg->channels == 0 || cfg->channels > FU_MAX_CH || x->k_max >= NBIN - LOBE - 1 ||
        x->k_top < 2 * k_floor || x->k_max <= x->k_harm) {
        free(x);
        return NULL;
    }
    x->fold_gain = pow(10.0, cfg->fold_db / 20.0);
    x->harm_gain = pow(10.0, cfg->harm_db / 20.0);
    x->skip = N - H;
    x->cap = N + 8192;
    x->frame = fftw_malloc(N * sizeof *x->frame);
    x->spec = fftw_malloc(NBIN * sizeof *x->spec);
    x->synth = fftw_malloc(NBIN * sizeof *x->synth);
    int ok = x->frame && x->spec && x->synth;
    for (uint32_t c = 0; c < cfg->channels && ok; c++) {
        x->ch[c].in = calloc(x->cap, sizeof(double));
        x->ch[c].ola = calloc(N, sizeof(double));
        x->ch[c].fill = N - H;
        ok = x->ch[c].in && x->ch[c].ola;
    }
    if (!ok) {
        fu_extend_free(x);
        return NULL;
    }
    x->fwd = fftw_plan_dft_r2c_1d(N, x->frame, x->spec, FFTW_ESTIMATE);
    x->inv = fftw_plan_dft_c2r_1d(N, x->synth, x->frame, FFTW_ESTIMATE);
    for (int i = 0; i < N; i++)
        x->win[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / N);
    return x;
}

void fu_extend_free(fu_extend *x)
{
    if (!x)
        return;
    if (x->fwd)
        fftw_destroy_plan(x->fwd);
    if (x->inv)
        fftw_destroy_plan(x->inv);
    fftw_free(x->frame);
    fftw_free(x->spec);
    fftw_free(x->synth);
    for (int c = 0; c < FU_MAX_CH; c++) {
        free(x->ch[c].in);
        free(x->ch[c].ola);
    }
    free(x);
}

const fu_extend_stats *fu_extend_stats_get(const fu_extend *x)
{
    return &x->st;
}

size_t fu_extend_bound(const fu_extend *x, size_t n_in)
{
    (void)x;
    return n_in + N + H;
}

/* Mean power over bins [a, b). */
static double band_mean(const double *cs, int a, int b)
{
    return (cs[b] - cs[a]) / (double)(b - a);
}

/* Folding: tile the normalised source band over [k_fold, k_fold_end). */
static void fold(fu_extend *x, chan *c)
{
    const double per_khz = 1000.0 / x->hz;
    const int kb = (int)per_khz;
    double lev[4], mean_l = 0.0;
    for (int j = 0; j < 4; j++) {
        const int b = x->k_top - (int)((4 - j) * per_khz);
        lev[j] = 10.0 * log10(band_mean(x->cs, b, b + kb) + 1e-300);
        mean_l += 0.25 * lev[j];
    }
    /* Least-squares slope over band centres -3.5, -2.5, -1.5, -0.5 kHz. */
    double slope = (-1.5 * (lev[0] - mean_l) - 0.5 * (lev[1] - mean_l) + 0.5 * (lev[2] - mean_l)
                    + 1.5 * (lev[3] - mean_l)) / 5.0;
    slope = fmin(fmax(slope, -12.0), -1.0);
    c->slope = c->primed ? 0.8 * c->slope + 0.2 * slope : slope;

    const double ref_db = lev[3];
    for (int k = x->k_fold; k < x->k_fold_end; k++) {
        const int ks = x->k_src + (k - x->k_fold) % x->src_bins;
        const int a = ks - SMOOTH < 0 ? 0 : ks - SMOOTH;
        const int b = ks + SMOOTH + 1 > x->k_top ? x->k_top : ks + SMOOTH + 1;
        const double dens = band_mean(x->cs, a, b);
        const double tgt_db = ref_db + c->slope * ((k - x->k_top) / per_khz + 0.5);
        double g = sqrt(pow(10.0, tgt_db / 10.0) / (dens + 1e-300));
        g = fmin(g, 4.0) * x->fold_gain;
        x->synth[k][0] += g * x->spec[ks][0];
        x->synth[k][1] += g * x->spec[ks][1];
    }
}

/* Harmonics: integer multiples of every tonal peak, placed above k_harm. */
static void harmonics(fu_extend *x, chan *c)
{
    const int half = x->k_top / 2, quarter = x->k_top / 4;
    const double lo = band_mean(x->cs, quarter, half), hi = band_mean(x->cs, half, x->k_top);
    c->dens_lo = c->primed ? 0.95 * c->dens_lo + 0.05 * lo : lo;
    c->dens_hi = c->primed ? 0.95 * c->dens_hi + 0.05 * hi : hi;
    const double tilt = fmin(fmax(10.0 * log10((c->dens_lo + 1e-300) / (c->dens_hi + 1e-300)), 6.0),
                             24.0);

    for (int k = x->k_peak0 > FLOOR_HALF ? x->k_peak0 : FLOOR_HALF; k < x->k_top - FLOOR_SKIP - 1; k++) {
        const double p = x->pw[k];
        if (!(p > x->pw[k - 1] && p >= x->pw[k + 1] && p >= x->pw[k - 2] && p >= x->pw[k + 2]))
            continue;
        const int a = k - FLOOR_HALF, b = k + FLOOR_HALF + 1 > x->k_top ? x->k_top : k + FLOOR_HALF + 1;
        const double around = (x->cs[b] - x->cs[a]) - (x->cs[k + FLOOR_SKIP + 1] - x->cs[k - FLOOR_SKIP]);
        const double floor_p = around / (double)(b - a - (2 * FLOOR_SKIP + 1));
        if (!(p > TONAL_RATIO * floor_p) || p < 1e-30)
            continue;

        /* Parabolic interpolation of the log magnitude: fractional bin d. */
        const double lm = 0.5 * log(x->pw[k - 1] + 1e-300), l0 = 0.5 * log(p),
                     lp = 0.5 * log(x->pw[k + 1] + 1e-300);
        const double den = lm - 2.0 * l0 + lp;
        const double d = den < 0.0 ? fmin(fmax(0.5 * (lm - lp) / den, -0.5), 0.5) : 0.0;
        const double kd = kernel(d);
        const double mag = sqrt(p) / kd;
        const double phi = atan2(x->spec[k][1], x->spec[k][0]);

        int wrote = 0;
        for (int h = 2; h <= x->cfg.max_harm; h++) {
            const double kappa = h * (k + d);
            if (kappa < x->k_harm + LOBE)
                continue;
            if (kappa > x->k_max)
                break;
            const double amp = mag * pow(10.0, -tilt * log2((double)h) / 20.0) * x->harm_gain;
            const double re = amp * cos(h * phi), im = amp * sin(h * phi);
            /* 2 * LOBE bins: the kernel's skirts are written too, or their absence
             * leaves spurs beside the partial (measured: 48 dB down with 4 bins). */
            const int j0 = (int)floor(kappa) - LOBE + 1;
            for (int j = j0; j < j0 + 2 * LOBE; j++) {
                const double w = kernel(j - kappa);
                x->synth[j][0] += w * re;
                x->synth[j][1] += w * im;
            }
            wrote++;
        }
        x->st.peaks += wrote > 0;
        x->st.lobes += (uint64_t)wrote;
    }
}

/* Analyse the N frames at c->in[off], overlap-add their extension into c->ola. */
static void process(fu_extend *x, chan *c, size_t off)
{
    for (int i = 0; i < N; i++)
        x->frame[i] = c->in[off + (size_t)i] * x->win[i];
    fftw_execute(x->fwd);

    /* Centre the phase reference; power and its prefix sums over the source band. */
    x->cs[0] = 0.0;
    for (int k = 0; k <= x->k_top; k++) {
        const double s = (k & 1) ? -1.0 : 1.0;
        x->spec[k][0] *= s;
        x->spec[k][1] *= s;
        x->pw[k] = x->spec[k][0] * x->spec[k][0] + x->spec[k][1] * x->spec[k][1];
        x->cs[k + 1] = x->cs[k] + x->pw[k];
    }
    x->st.frames++;
    if (x->cs[x->k_top] < 1e-24)
        return;

    memset(x->synth, 0, NBIN * sizeof *x->synth);
    if (x->src_bins >= 8)
        fold(x, c);
    harmonics(x, c);
    c->primed = 1;

    /* Tweeter guard: what is synthesised above fold_top may not exceed
     * ULTRA_CAP of the real signal's power in this frame. */
    double ultra = 0.0;
    for (int k = x->k_harm; k <= x->k_max + LOBE; k++)
        ultra += x->synth[k][0] * x->synth[k][0] + x->synth[k][1] * x->synth[k][1];
    const double allowed = ULTRA_CAP * x->cs[x->k_top];
    const double trim = ultra > allowed ? sqrt(allowed / ultra) : 1.0;
    x->st.capped += trim < 1.0;

    for (int k = x->k_fold; k <= x->k_max + LOBE; k++) {
        const double s = ((k & 1) ? -1.0 : 1.0) * (k >= x->k_harm ? trim : 1.0);
        x->synth[k][0] *= s;
        x->synth[k][1] *= s;
    }
    fftw_execute(x->inv);
    for (int i = 0; i < N; i++)
        c->ola[i] += x->frame[i] * x->win[i] * (OLA_NORM / N);
}

/* Append n frames (zeros when in == NULL); emit every completed hop. */
static size_t push(fu_extend *x, const double *const in[], size_t n, double *const out[],
                   size_t out_cap, uint64_t limit)
{
    const uint32_t nch = x->cfg.channels;
    if (x->ch[0].fill + n > x->cap) {
        const size_t cap = 2 * (x->ch[0].fill + n);
        for (uint32_t c = 0; c < nch; c++) {
            double *b = realloc(x->ch[c].in, cap * sizeof *b);
            if (!b)
                return 0;
            x->ch[c].in = b;
        }
        x->cap = cap;
    }
    size_t emitted = 0;
    for (uint32_t c = 0; c < nch; c++) {
        chan *ch = &x->ch[c];
        if (in)
            memcpy(ch->in + ch->fill, in[c], n * sizeof(double));
        else
            memset(ch->in + ch->fill, 0, n * sizeof(double));
        ch->fill += n;

        size_t off = 0, w = 0, skip = x->skip;
        uint64_t done = x->total_out;
        while (ch->fill - off >= N) {
            process(x, ch, off);
            for (int i = 0; i < H; i++) {
                if (skip) {
                    skip--;
                } else if (done < limit && w < out_cap) {
                    out[c][w++] = ch->in[off + (size_t)i] + ch->ola[i];
                    done++;
                }
            }
            memmove(ch->ola, ch->ola + H, (N - H) * sizeof(double));
            memset(ch->ola + N - H, 0, H * sizeof(double));
            off += H;
        }
        memmove(ch->in, ch->in + off, (ch->fill - off) * sizeof(double));
        ch->fill -= off;
        if (c == nch - 1) {
            x->skip = skip;
            emitted = w;
        }
    }
    x->total_out += emitted;
    return emitted;
}

size_t fu_extend_run(fu_extend *x, const double *const in[], size_t n_in,
                     double *const out[], size_t out_cap)
{
    x->total_in += n_in;
    return push(x, in, n_in, out, out_cap, UINT64_MAX);
}

size_t fu_extend_flush(fu_extend *x, double *const out[], size_t out_cap)
{
    /* Enough zeros to complete the last hop and drain the N - H frames of delay. */
    const size_t pad = (N - H) + (H - (size_t)(x->total_in % H)) % H;
    return push(x, NULL, pad, out, out_cap, x->total_in);
}
