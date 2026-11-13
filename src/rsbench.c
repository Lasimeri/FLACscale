/*
 * rsbench.c -- gate S1: measure the resampler before anything is built on it.
 *
 * Purpose : four tables, all from this machine, all from synthetic signals
 *           (white noise and pure tones, generated here, labelled as such):
 *             1. filter response of each stage design (passband ripple,
 *                stopband peak), from the DTFT of the prototype taps
 *             2. correctness: production kernels against the scalar formula,
 *                and chunked streaming against itself
 *             3. tone test: gain, residual and time offset at 192 kHz,
 *                FLACscale against libsoxr
 *             4. speed: 3 interleaved rounds of every variant
 * Usage   : rsbench            (no arguments; about 20 s)
 * Baseline: libsoxr 0.1.3, recipes SOXR_VHQ and SOXR_32_BITQ, 1 thread.
 */
#include "fu_resample.h"
#include "fu_resample_ref.h"

#include <math.h>
#include <soxr.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum { RATE_IN = 44100, RATE_OUT = 192000, L1 = 640, M1 = 147, NSTEEP = 256, NRELAX = 32 };
#define PASS 0.907

static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return 1e3 * (double)ts.tv_sec + 1e-6 * (double)ts.tv_nsec;
}

static void *xalloc(size_t n, size_t size)
{
    void *p = calloc(n, size);
    if (!p) {
        fprintf(stderr, "rsbench: out of memory (%zu x %zu bytes)\n", n, size);
        exit(1);
    }
    return p;
}

/* ---- 1. filter response ------------------------------------------------ */

/* |H| in dB at f cycles per input sample; h runs at L times the input rate. */
static double mag_db(const double *h, size_t len, uint32_t L, double f)
{
    const long double w = 2.0L * (long double)M_PI * (long double)f / (long double)L;
    const double dr = (double)cosl(w), di = (double)sinl(w);
    double re = 0.0, im = 0.0;
    for (size_t k0 = 0; k0 < len; k0 += 256) {
        double cr = (double)cosl(w * (long double)k0), ci = (double)sinl(w * (long double)k0);
        const size_t k1 = k0 + 256 < len ? k0 + 256 : len;
        for (size_t k = k0; k < k1; k++) {
            re += h[k] * cr;
            im -= h[k] * ci;
            const double t = cr * dr - ci * di;
            ci = cr * di + ci * dr;
            cr = t;
        }
    }
    return 20.0 * log10(hypot(re, im) / (double)L);
}

static void response_row(const char *name, uint32_t L, uint32_t N, double fpass, double fstop,
                         double rate_in)
{
    const size_t len = (size_t)L * N;
    double *h = xalloc(len, sizeof *h);
    fu_design(L, N, fpass, fstop, h);

    double lo = 1e9, hi = -1e9;
    for (int i = 0; i <= 1000; i++) {
        const double d = mag_db(h, len, L, fpass * i / 1000.0);
        lo = d < lo ? d : lo;
        hi = d > hi ? d : hi;
    }
    /* Kaiser side lobes fall away from the edge: scan the first 40 lobes
     * densely, then the rest of the band up to the high-rate Nyquist. */
    double stop = -1e9;
    const double dense = 40.0 / N;
    for (int i = 0; i <= 320; i++) {
        const double d = mag_db(h, len, L, fstop + dense * i / 320.0);
        stop = d > stop ? d : stop;
    }
    const double f0 = fstop + dense, f1 = 0.5 * L;
    for (int i = 0; i <= 2000 && f1 > f0; i++) {
        const double d = mag_db(h, len, L, f0 * pow(f1 / f0, i / 2000.0));
        stop = d > stop ? d : stop;
    }
    printf("  %-22s %4u %4u %8.0f %8.0f %9.1f %12.3e %10.1f\n", name, L, N,
           fpass * rate_in, fstop * rate_in, fu_design_atten(N, fpass, fstop), hi - lo, stop);
    free(h);
}

/* ---- runners ----------------------------------------------------------- */

static size_t run_stage(uint32_t L, uint32_t M, uint32_t N, double fpass, double fstop,
                        const double *xl, const double *xr, size_t n_in,
                        double *yl, double *yr, size_t cap, const size_t *chunks, size_t nchunks)
{
    fu_stage s;
    if (fu_stage_init(&s, L, M, N, fpass, fstop, 2)) {
        fprintf(stderr, "rsbench: stage init failed\n");
        exit(1);
    }
    size_t done = 0, n = 0;
    for (size_t i = 0; done < n_in; i++) {
        size_t take = chunks[i % nchunks];
        take = take < n_in - done ? take : n_in - done;
        const double *in[2] = { xl + done, xr + done };
        double *out[2] = { yl + n, yr + n };
        n += fu_stage_run(&s, in, take, out, cap - n);
        done += take;
    }
    double *out[2] = { yl + n, yr + n };
    n += fu_stage_flush(&s, out, cap - n);
    fu_stage_free(&s);
    return n;
}

static size_t run_chain(uint32_t nch, const double *xl, const double *xr, size_t n_in,
                        double *yl, double *yr, size_t cap)
{
    fu_chain c;
    if (fu_chain_init(&c, RATE_IN, RATE_OUT, nch, NULL)) {
        fprintf(stderr, "rsbench: chain init failed\n");
        exit(1);
    }
    size_t done = 0, n = 0;
    while (done < n_in) {
        const size_t take = n_in - done < 4096 ? n_in - done : 4096;
        const double *in[2] = { xl + done, xr ? xr + done : NULL };
        double *out[2] = { yl + n, yr ? yr + n : NULL };
        n += fu_chain_run(&c, in, take, out, cap - n);
        done += take;
    }
    double *out[2] = { yl + n, yr ? yr + n : NULL };
    n += fu_chain_flush(&c, out, cap - n);
    fu_chain_free(&c);
    return n;
}

static size_t run_soxr(unsigned long recipe, uint32_t nch, const double *xl, const double *xr,
                       size_t n_in, double *yl, double *yr, size_t cap)
{
    const soxr_quality_spec_t q  = soxr_quality_spec(recipe, 0);
    const soxr_io_spec_t      io = soxr_io_spec(SOXR_FLOAT64_S, SOXR_FLOAT64_S);
    const soxr_runtime_spec_t rt = soxr_runtime_spec(1);
    const double *in[2] = { xl, xr };
    double *out[2] = { yl, yr };
    size_t idone = 0, odone = 0;
    const soxr_error_t e = soxr_oneshot(RATE_IN, RATE_OUT, nch, in, n_in, &idone,
                                        out, cap, &odone, &io, &q, &rt);
    if (e) {
        fprintf(stderr, "rsbench: soxr: %s\n", soxr_strerror(e));
        exit(1);
    }
    return odone;
}

static double max_abs_diff(const double *a, const double *b, size_t n)
{
    double m = 0.0;
    for (size_t i = 0; i < n; i++) {
        const double d = fabs(a[i] - b[i]);
        m = d > m ? d : m;
    }
    return m;
}

/* ---- 3. tone test ------------------------------------------------------ */

typedef struct tone_fit {
    double gain_db, resid_db, offset_ns;
} tone_fit;

/* Least-squares fit of a*sin + b*cos at f over y[n0..n1) sampled at fs. */
static tone_fit fit_tone(const double *y, size_t n0, size_t n1, double f, double fs, double amp)
{
    const double w = 2.0 * M_PI * f / fs;
    double ss = 0, sc = 0, cc = 0, ys = 0, yc = 0;
    for (size_t n = n0; n < n1; n++) {
        const double s = sin(w * (double)n), c = cos(w * (double)n);
        ss += s * s;
        sc += s * c;
        cc += c * c;
        ys += y[n] * s;
        yc += y[n] * c;
    }
    const double det = ss * cc - sc * sc;
    const double a = (ys * cc - yc * sc) / det, b = (yc * ss - ys * sc) / det;
    double er = 0, ef = 0;
    for (size_t n = n0; n < n1; n++) {
        const double fit = a * sin(w * (double)n) + b * cos(w * (double)n);
        er += (y[n] - fit) * (y[n] - fit);
        ef += fit * fit;
    }
    return (tone_fit){
        .gain_db   = 20.0 * log10(hypot(a, b) / amp),
        .resid_db  = 10.0 * log10(er / ef + 1e-40),
        .offset_ns = -atan2(b, a) / (2.0 * M_PI * f) * 1e9,
    };
}

static void tone_table(void)
{
    static const double freqs[] = { 30.0, 1000.0, 10000.0, 19000.0, 20000.0 };
    const size_t n_in = RATE_IN, cap = RATE_OUT + 4096;
    double *x = xalloc(n_in, sizeof *x), *y = xalloc(cap, sizeof *y);

    printf("\n3. Tone test: synthetic sine, amplitude 0.5, 1 s, fitted over 0.1 s .. 0.9 s at 192 kHz\n");
    printf("  %8s | %-36s | %-36s\n", "", "FLACscale two-stage", "libsoxr VHQ");
    printf("  %8s | %11s %11s %11s | %11s %11s %11s\n", "tone Hz", "gain dB", "resid dB",
           "offset ns", "gain dB", "resid dB", "offset ns");
    for (size_t i = 0; i < sizeof freqs / sizeof *freqs; i++) {
        for (size_t k = 0; k < n_in; k++)
            x[k] = 0.5 * sin(2.0 * M_PI * freqs[i] * (double)k / RATE_IN);
        const size_t n0 = RATE_OUT / 10, n1 = 9 * (RATE_OUT / 10);
        run_chain(1, x, NULL, n_in, y, NULL, cap);
        const tone_fit a = fit_tone(y, n0, n1, freqs[i], RATE_OUT, 0.5);
        run_soxr(SOXR_VHQ, 1, x, NULL, n_in, y, NULL, cap);
        const tone_fit b = fit_tone(y, n0, n1, freqs[i], RATE_OUT, 0.5);
        printf("  %8.0f | %11.2e %11.1f %11.3f | %11.2e %11.1f %11.3f\n", freqs[i],
               a.gain_db, a.resid_db, a.offset_ns, b.gain_db, b.resid_db, b.offset_ns);
    }
    free(x);
    free(y);
}

/* ---- main -------------------------------------------------------------- */

int main(void)
{
    const double fp = 0.5 * PASS;

    printf("1. Filter response (DTFT of the prototype taps; frequencies in Hz at the stage input)\n");
    printf("  %-22s %4s %4s %8s %8s %9s %12s %10s\n", "design", "L", "N", "pass Hz", "stop Hz",
           "design dB", "ripple dB pp", "stop dB");
    response_row("single 640/147", L1, NSTEEP, fp, 0.5, RATE_IN);
    response_row("stage A 2/1 steep", 2, NSTEEP, fp, 0.5, RATE_IN);
    response_row("stage B 320/147 relax", 320, NRELAX, 0.25, 0.75, 2.0 * RATE_IN);

    /* Synthetic input: 10 s of stereo white noise, uniform in [-0.5, 0.5). */
    const size_t n_in = 10 * (size_t)RATE_IN;
    const size_t n_out = fu_ref_count(n_in, L1, M1), cap = n_out + 4096;
    double *xl = xalloc(n_in, sizeof *xl), *xr = xalloc(n_in, sizeof *xr);
    double *ol = xalloc(cap, sizeof *ol), *or_ = xalloc(cap, sizeof *or_);
    double *yl = xalloc(cap, sizeof *yl), *yr = xalloc(cap, sizeof *yr);
    double *zl = xalloc(cap, sizeof *zl), *zr = xalloc(cap, sizeof *zr);
    uint64_t rng = 0x9E3779B97F4A7C15ull;
    for (size_t k = 0; k < 2 * n_in; k++) {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        (k & 1 ? xr : xl)[k >> 1] = (double)(rng >> 11) * 0x1p-53 - 0.5;
    }
    double *h = xalloc((size_t)L1 * NSTEEP, sizeof *h);
    fu_design(L1, NSTEEP, fp, 0.5, h);

    static const size_t even[] = { 4096 };
    static const size_t odd[]  = { 1, 4096, 7, 777, 2, 65536, 255, 3 };

    printf("\n2. Correctness on 10 s synthetic stereo noise (%zu frames in, %zu expected out)\n",
           n_in, n_out);
    const size_t c0 = fu_ref_scalar(L1, M1, NSTEEP, h, xl, n_in, ol);
    fu_ref_scalar(L1, M1, NSTEEP, h, xr, n_in, or_);
    const size_t c1 = fu_ref_phase2(L1, M1, NSTEEP, h, xl, xr, n_in, yl, yr);
    printf("  %-44s frames %zu  max|d| vs scalar formula %.3e\n", "AVX2 phase-order (plain layout)",
           c1, fmax(max_abs_diff(ol, yl, c0), max_abs_diff(or_, yr, c0)));
    const size_t c2 = run_stage(L1, M1, NSTEEP, fp, 0.5, xl, xr, n_in, yl, yr, cap, even, 1);
    printf("  %-44s frames %zu  max|d| vs scalar formula %.3e\n", "AVX2 visit-order stage, 4096-frame chunks",
           c2, fmax(max_abs_diff(ol, yl, c0), max_abs_diff(or_, yr, c0)));
    const size_t c3 = run_stage(L1, M1, NSTEEP, fp, 0.5, xl, xr, n_in, zl, zr, cap, odd, 8);
    printf("  %-44s frames %zu  max|d| vs 4096-frame chunks %.3e\n", "AVX2 visit-order stage, ragged chunks",
           c3, fmax(max_abs_diff(yl, zl, c2), max_abs_diff(yr, zr, c2)));
    const size_t c4 = run_chain(2, xl, xr, n_in, zl, zr, cap);
    printf("  %-44s frames %zu  (different filter: compared in table 3)\n", "two-stage chain", c4);
    if (c0 != n_out || c1 != n_out || c2 != n_out || c3 != n_out || c4 != n_out) {
        printf("  FAIL: frame count differs from ceil(n * 640 / 147) = %zu\n", n_out);
        return 1;
    }

    tone_table();

    printf("\n4. Speed, 10 s stereo, 44100 -> 192000, 1 thread, 3 interleaved rounds\n");
    enum { NV = 6, ROUNDS = 3 };
    static const char *const name[NV] = {
        "scalar formula, single stage 640/147 x256",
        "AVX2 phase-order, single stage",
        "AVX2 visit-order, single stage",
        "AVX2 two-stage chain (2/1 x256, 320/147 x32)",
        "libsoxr VHQ (28-bit recipe)",
        "libsoxr 32-bit recipe",
    };
    double best[NV], worst[NV];
    for (int v = 0; v < NV; v++) {
        best[v] = 1e30;
        worst[v] = 0.0;
    }
    for (int r = 0; r < ROUNDS; r++) {
        for (int v = 0; v < NV; v++) {
            const double t0 = now_ms();
            switch (v) {
            case 0:
                fu_ref_scalar(L1, M1, NSTEEP, h, xl, n_in, yl);
                fu_ref_scalar(L1, M1, NSTEEP, h, xr, n_in, yr);
                break;
            case 1: fu_ref_phase2(L1, M1, NSTEEP, h, xl, xr, n_in, yl, yr); break;
            case 2: run_stage(L1, M1, NSTEEP, fp, 0.5, xl, xr, n_in, yl, yr, cap, even, 1); break;
            case 3: run_chain(2, xl, xr, n_in, yl, yr, cap); break;
            case 4: run_soxr(SOXR_VHQ, 2, xl, xr, n_in, yl, yr, cap); break;
            default: run_soxr(SOXR_32_BITQ, 2, xl, xr, n_in, yl, yr, cap); break;
            }
            const double ms = now_ms() - t0;
            best[v] = ms < best[v] ? ms : best[v];
            worst[v] = ms > worst[v] ? ms : worst[v];
        }
    }
    printf("  %-46s %10s %10s %12s\n", "variant", "best ms", "worst ms", "x real time");
    for (int v = 0; v < NV; v++)
        printf("  %-46s %10.1f %10.1f %12.1f\n", name[v], best[v], worst[v], 10000.0 / best[v]);

    free(h);
    free(xl); free(xr); free(ol); free(or_);
    free(yl); free(yr); free(zl); free(zr);
    return 0;
}
