/*
 * fu_dot.h -- the one inner product every resampler stage runs.
 *
 * Purpose : sum(c[i] * x[i]) over N taps in double precision, for one channel
 *           (fu_dot1) or for two channels against the same coefficients
 *           (fu_dot2, one coefficient load feeds both channels).
 * Contract: N is a multiple of 16; c is 32-byte aligned; x may be unaligned.
 *           Four accumulators per channel: on Zen 3 an FMA has latency 4 and
 *           two ports, so eight independent chains (2 channels x 4) keep both
 *           ports busy every cycle.
 * Target  : AVX2 + FMA3 (Ryzen 7 5800X).
 */
#ifndef FU_DOT_H
#define FU_DOT_H

#include <immintrin.h>
#include <stdint.h>

static inline double fu_hsum(__m256d v)
{
    __m128d lo = _mm_add_pd(_mm256_castpd256_pd128(v), _mm256_extractf128_pd(v, 1));
    return _mm_cvtsd_f64(_mm_add_sd(lo, _mm_unpackhi_pd(lo, lo)));
}

[[maybe_unused]] static inline double fu_dot1(const double *restrict c, const double *restrict x, uint32_t N)
{
    __m256d a0 = _mm256_setzero_pd(), a1 = a0, a2 = a0, a3 = a0;
    for (uint32_t i = 0; i < N; i += 16) {
        a0 = _mm256_fmadd_pd(_mm256_load_pd(c + i),      _mm256_loadu_pd(x + i),      a0);
        a1 = _mm256_fmadd_pd(_mm256_load_pd(c + i + 4),  _mm256_loadu_pd(x + i + 4),  a1);
        a2 = _mm256_fmadd_pd(_mm256_load_pd(c + i + 8),  _mm256_loadu_pd(x + i + 8),  a2);
        a3 = _mm256_fmadd_pd(_mm256_load_pd(c + i + 12), _mm256_loadu_pd(x + i + 12), a3);
    }
    return fu_hsum(_mm256_add_pd(_mm256_add_pd(a0, a1), _mm256_add_pd(a2, a3)));
}

[[maybe_unused]] static inline void fu_dot2(const double *restrict c, const double *restrict xl,
                           const double *restrict xr, uint32_t N,
                           double *restrict yl, double *restrict yr)
{
    __m256d l0 = _mm256_setzero_pd(), l1 = l0, l2 = l0, l3 = l0;
    __m256d r0 = l0, r1 = l0, r2 = l0, r3 = l0;
    for (uint32_t i = 0; i < N; i += 16) {
        const __m256d c0 = _mm256_load_pd(c + i);
        const __m256d c1 = _mm256_load_pd(c + i + 4);
        const __m256d c2 = _mm256_load_pd(c + i + 8);
        const __m256d c3 = _mm256_load_pd(c + i + 12);
        l0 = _mm256_fmadd_pd(c0, _mm256_loadu_pd(xl + i),      l0);
        r0 = _mm256_fmadd_pd(c0, _mm256_loadu_pd(xr + i),      r0);
        l1 = _mm256_fmadd_pd(c1, _mm256_loadu_pd(xl + i + 4),  l1);
        r1 = _mm256_fmadd_pd(c1, _mm256_loadu_pd(xr + i + 4),  r1);
        l2 = _mm256_fmadd_pd(c2, _mm256_loadu_pd(xl + i + 8),  l2);
        r2 = _mm256_fmadd_pd(c2, _mm256_loadu_pd(xr + i + 8),  r2);
        l3 = _mm256_fmadd_pd(c3, _mm256_loadu_pd(xl + i + 12), l3);
        r3 = _mm256_fmadd_pd(c3, _mm256_loadu_pd(xr + i + 12), r3);
    }
    *yl = fu_hsum(_mm256_add_pd(_mm256_add_pd(l0, l1), _mm256_add_pd(l2, l3)));
    *yr = fu_hsum(_mm256_add_pd(_mm256_add_pd(r0, r1), _mm256_add_pd(r2, r3)));
}

#endif
