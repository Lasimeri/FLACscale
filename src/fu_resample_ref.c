/*
 * fu_resample_ref.c -- plain reference forms (see fu_resample_ref.h).
 *
 * Math    : y[n] = sum_j h[p + j*L] * x[q - j],  n*M + L*N/2 = q*L + p,
 *           x[k] = 0 outside 0 <= k < n_in.
 */
#include "fu_resample_ref.h"
#include "fu_dot.h"

#include <stdlib.h>
#include <string.h>

size_t fu_ref_count(size_t n_in, uint32_t L, uint32_t M)
{
    return (size_t)(((uint64_t)n_in * L + M - 1) / M);
}

size_t fu_ref_scalar(uint32_t L, uint32_t M, uint32_t N, const double *h,
                     const double *x, size_t n_in, double *y)
{
    const size_t n_out = fu_ref_count(n_in, L, M);
    const uint64_t c = (uint64_t)L * N / 2;

    for (size_t n = 0; n < n_out; n++) {
        const uint64_t t = (uint64_t)n * M + c;
        const int64_t  q = (int64_t)(t / L);
        const uint32_t p = (uint32_t)(t % L);
        double acc = 0.0;
        for (uint32_t j = 0; j < N; j++) {
            const int64_t k = q - (int64_t)j;
            if (k >= 0 && k < (int64_t)n_in)
                acc += h[p + (size_t)j * L] * x[k];
        }
        y[n] = acc;
    }
    return n_out;
}

size_t fu_ref_phase2(uint32_t L, uint32_t M, uint32_t N, const double *h,
                     const double *xl, const double *xr, size_t n_in,
                     double *yl, double *yr)
{
    const size_t n_out = fu_ref_count(n_in, L, M);
    const size_t pre = N / 2 - 1, post = N / 2;
    const size_t len = pre + n_in + post;
    double *tab = aligned_alloc(64, ((size_t)L * N * sizeof *tab + 63) & ~(size_t)63);
    double *bl = calloc(len, sizeof *bl);
    double *br = calloc(len, sizeof *br);
    if (!tab || !bl || !br) {
        free(tab);
        free(bl);
        free(br);
        return 0;
    }
    for (uint32_t p = 0; p < L; p++)
        for (uint32_t i = 0; i < N; i++)
            tab[(size_t)p * N + i] = h[p + (size_t)(N - 1 - i) * L];
    memcpy(bl + pre, xl, n_in * sizeof *bl);
    memcpy(br + pre, xr, n_in * sizeof *br);

    const uint64_t c = (uint64_t)L * N / 2;
    for (size_t n = 0; n < n_out; n++) {
        const uint64_t t = (uint64_t)n * M + c;
        const size_t   w = (size_t)(t / L) - N / 2;
        const uint32_t p = (uint32_t)(t % L);
        fu_dot2(tab + (size_t)p * N, bl + w, br + w, N, yl + n, yr + n);
    }
    free(tab);
    free(bl);
    free(br);
    return n_out;
}
