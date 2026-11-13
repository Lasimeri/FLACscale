/*
 * fu_dequant.c -- consistent dequantization (see fu_dequant.h).
 *
 * Frame    : periodic Hann of N = 2048 scaled by 1/sqrt(1.5), hop H = N/4, with
 *            a unitary FFT. The squares of the shifted windows sum to 1, so
 *            synthesis after analysis is the identity. Frames wrap around the
 *            segment, which makes that exact at the segment edges; the seam
 *            lies in the margin, which is thrown away.
 * Segments : SEG samples are solved at a time; the first and last MARGIN of
 *            them are context. Neighbouring results are blended linearly over
 *            2 * FADE samples. A blend of two consistent values is consistent.
 * Noise    : a coefficient of this transform carries a quarter of the
 *            per-sample rounding-noise power (the window: 0.375 / 1.5).
 */
#include "fu_dequant.h"

#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

enum {
    N = 2048, H = N / 4, NBIN = N / 2 + 1,
    SEG = 1 << 16, FRAMES = SEG / H,
    MARGIN = 2 * N, FADE = N / 4,
    STRIDE = SEG - 2 * MARGIN, SPAN = STRIDE + 2 * FADE
};

typedef struct workspace {
    fftw_complex *c;        /* FRAMES * NBIN coefficients of the rounded signal */
    fftw_complex *spec;     /* one frame's spectrum */
    double       *frame;    /* one frame in time */
    double       *q;        /* the rounded segment */
    double       *u;        /* the estimate */
    double       *lo, *hi;  /* the interval of every sample */
    float        *p, *g;    /* per coefficient: power, gain */
} workspace;

typedef struct transform {
    const double *win;
    fftw_plan     fwd, inv;
    double        unit;     /* 1/sqrt(N): makes each FFT unitary */
} transform;

fu_dequant_cfg fu_dequant_cfg_default(void)
{
    return (fu_dequant_cfg){ .passes = 2, .half = 0.5, .noise = 1.0 / 12.0, .threads = 0 };
}

static void workspace_free(workspace *w)
{
    fftw_free(w->c);
    fftw_free(w->spec);
    fftw_free(w->frame);
    fftw_free(w->q);
    fftw_free(w->u);
    fftw_free(w->lo);
    fftw_free(w->hi);
    fftw_free(w->p);
    fftw_free(w->g);
}

static int workspace_init(workspace *w)
{
    w->c = fftw_alloc_complex((size_t)FRAMES * NBIN);
    w->spec = fftw_alloc_complex(NBIN);
    w->frame = fftw_alloc_real(N);
    w->q = fftw_alloc_real(SEG);
    w->u = fftw_alloc_real(SEG);
    w->lo = fftw_alloc_real(SEG);
    w->hi = fftw_alloc_real(SEG);
    w->p = fftw_malloc((size_t)FRAMES * NBIN * sizeof(float));
    w->g = fftw_malloc((size_t)FRAMES * NBIN * sizeof(float));
    return w->c && w->spec && w->frame && w->q && w->u && w->lo && w->hi && w->p && w->g ? 0 : -1;
}

/* Spectrum of frame k of the segment signal v, into w->spec (not yet scaled). */
static void analyse(const workspace *w, const transform *t, const double *v, int k)
{
    const int at = k * H;
    for (int j = 0; j < N; j++)
        w->frame[j] = t->win[j] * v[(at + j) & (SEG - 1)];
    fftw_execute_dft_r2c(t->fwd, w->frame, w->spec);
}

/* w->p = power of every coefficient of v; keep the coefficients in w->c when asked. */
static void powers(const workspace *w, const transform *t, const double *v, int keep)
{
    for (int k = 0; k < FRAMES; k++) {
        analyse(w, t, v, k);
        fftw_complex *c = w->c + (size_t)k * NBIN;
        float        *p = w->p + (size_t)k * NBIN;
        for (int b = 0; b < NBIN; b++) {
            const double re = t->unit * w->spec[b][0], im = t->unit * w->spec[b][1];
            p[b] = (float)(re * re + im * im);
            if (keep) {
                c[b][0] = re;
                c[b][1] = im;
            }
        }
    }
}

/* Wiener gain of every coefficient: S / (S + noise). S is the signal power:
 * w->p averaged over the 3 x 3 neighbouring frames and bins, less `known`
 * (the noise that the averaged power is known to contain). */
static void gains(const workspace *w, double noise, double known)
{
    for (int k = 0; k < FRAMES; k++)
        for (int b = 0; b < NBIN; b++) {
            double sum = 0.0;
            int    cells = 0;
            for (int dk = -1; dk <= 1; dk++) {
                const float *row = w->p + (size_t)((k + dk) & (FRAMES - 1)) * NBIN;
                for (int nb = b > 0 ? b - 1 : b; nb <= (b < NBIN - 1 ? b + 1 : b); nb++, cells++)
                    sum += row[nb];
            }
            const double s = fmax(sum / cells - known, 0.0);
            w->g[(size_t)k * NBIN + b] = (float)(s / (s + noise));
        }
}

/* w->u = the filtered rounded signal (gains applied to w->c), projected into the box. */
static void estimate(const workspace *w, const transform *t)
{
    memset(w->u, 0, SEG * sizeof(double));
    for (int k = 0; k < FRAMES; k++) {
        const fftw_complex *c = w->c + (size_t)k * NBIN;
        const float        *g = w->g + (size_t)k * NBIN;
        for (int b = 0; b < NBIN; b++) {
            w->spec[b][0] = g[b] * c[b][0];
            w->spec[b][1] = g[b] * c[b][1];
        }
        fftw_execute_dft_c2r(t->inv, w->spec, w->frame);
        const int at = k * H;
        for (int j = 0; j < N; j++)
            w->u[(at + j) & (SEG - 1)] += t->win[j] * w->frame[j];
    }
    for (int i = 0; i < SEG; i++)
        w->u[i] = fmin(fmax(t->unit * w->u[i], w->lo[i]), w->hi[i]);
}

/* One segment: Wiener filter of the rounded signal, then the projection. A
 * further pass re-estimates the signal power from the cleaner estimate and
 * filters the rounded signal again. */
static void solve(const workspace *w, const transform *t, double noise, int passes)
{
    powers(w, t, w->q, 1);
    gains(w, noise, noise);
    estimate(w, t);
    for (int pass = 1; pass < passes; pass++) {
        powers(w, t, w->u, 0);
        gains(w, noise, 0.0);
        estimate(w, t);
    }
}

int fu_dequant(double *x, size_t n, double step, const fu_dequant_cfg *cfg)
{
    if (n == 0)
        return 0;
    const size_t segments = (n + STRIDE - 1) / STRIDE;
    const int    threads = cfg->threads > 0 ? cfg->threads : (int)sysconf(_SC_NPROCESSORS_ONLN);
    double      *win = fftw_alloc_real(N);
    double      *res = malloc(segments * SPAN * sizeof *res);
    double      *probe_t = fftw_alloc_real(N);
    fftw_complex *probe_f = fftw_alloc_complex(NBIN);
    if (!win || !res || !probe_t || !probe_f) {
        fftw_free(win);
        free(res);
        fftw_free(probe_t);
        fftw_free(probe_f);
        return -1;
    }
    for (int j = 0; j < N; j++)
        win[j] = (0.5 - 0.5 * cos(2.0 * M_PI * j / N)) / sqrt(1.5);
    /* Plans are made once, serially (the planner is not thread-safe); executing
     * them on other arrays from several threads is. */
    const transform t = { .win = win, .unit = 1.0 / sqrt((double)N),
                          .fwd = fftw_plan_dft_r2c_1d(N, probe_t, probe_f, FFTW_MEASURE),
                          .inv = fftw_plan_dft_c2r_1d(N, probe_f, probe_t, FFTW_MEASURE) };
    const double half = cfg->half * step, noise = 0.25 * cfg->noise * step * step;
    int failed = 0;
#pragma omp parallel num_threads(threads)
    {
        workspace w = { 0 };
        const int bad = workspace_init(&w);
#pragma omp atomic
        failed += bad;
#pragma omp for schedule(dynamic)
        for (size_t s = 0; s < segments; s++) {
            if (bad)
                continue;
            /* Segment sample i is signal sample s * STRIDE - MARGIN + i. Outside
             * the signal the value is known to be zero. */
            const long long base = (long long)(s * STRIDE) - MARGIN;
            for (int i = 0; i < SEG; i++) {
                const long long at = base + i;
                const int inside = at >= 0 && (size_t)at < n;
                w.q[i] = inside ? x[at] : 0.0;
                w.lo[i] = inside ? w.q[i] - half : 0.0;
                w.hi[i] = inside ? w.q[i] + half : 0.0;
            }
            solve(&w, &t, noise, cfg->passes);
            memcpy(res + s * SPAN, w.u + MARGIN - FADE, SPAN * sizeof(double));
        }
        workspace_free(&w);
    }
    fftw_destroy_plan(t.fwd);
    fftw_destroy_plan(t.inv);
    fftw_free(probe_t);
    fftw_free(probe_f);
    fftw_free(win);
    if (failed) {
        free(res);
        return -1;
    }

    /* Blend. Result s covers signal samples [s*STRIDE - FADE, (s+1)*STRIDE + FADE);
     * its weight ramps 0 -> 1 over the first 2*FADE and 1 -> 0 over the last. */
    memset(x, 0, n * sizeof *x);
    for (size_t s = 0; s < segments; s++) {
        const long long first = (long long)(s * STRIDE) - FADE;
        const double *r = res + s * SPAN;
        for (int i = 0; i < SPAN; i++) {
            const long long at = first + i;
            if (at < 0 || (size_t)at >= n)
                continue;
            double wgt = 1.0;
            if (i < 2 * FADE && s > 0)
                wgt = (i + 0.5) / (2.0 * FADE);
            else if (i >= SPAN - 2 * FADE && s + 1 < segments)
                wgt = (SPAN - i - 0.5) / (2.0 * FADE);
            x[at] += wgt * r[i];
        }
    }
    free(res);
    return 0;
}
