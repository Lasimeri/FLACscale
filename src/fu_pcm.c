/*
 * fu_pcm.c -- double <-> integer PCM boundary (see fu_pcm.h).
 *
 * Technique: saturation and clip counting without a branch per sample. The
 *            clip test is two ordered compares whose sign masks are counted
 *            with popcount; the saturation is max/min in the double domain
 *            (both bounds are exact doubles), so vcvtpd2dq never sees a value
 *            it cannot represent.
 * Clip rule: with FS = 2^(bits-1): rounded value > FS - 1  <=>  v >= FS - 0.5
 *            (the tie rounds to the even FS); rounded value < -FS  <=>
 *            v < -FS - 0.5 (the tie rounds to the even -FS, which is in range).
 */
#include "fu_pcm.h"

#include <immintrin.h>
#include <math.h>

double fu_pcm_limit(unsigned bits)
{
    const double fs = (double)(1ull << (bits - 1));
    return (fs - 1.0) / fs;
}

size_t fu_pcm_to_int_ref(const double *x, int32_t *y, size_t n, double gain, unsigned bits)
{
    const double fs = (double)(1ull << (bits - 1)), k = gain * fs;
    size_t clipped = 0;
    for (size_t i = 0; i < n; i++) {
        const double v = x[i] * k;
        clipped += (v >= fs - 0.5) | (v < -fs - 0.5);
        y[i] = (int32_t)llrint(fmin(fmax(v, -fs), fs - 1.0));
    }
    return clipped;
}

size_t fu_pcm_to_int(const double *x, int32_t *y, size_t n, double gain, unsigned bits)
{
    const double  fs = (double)(1ull << (bits - 1));
    const __m256d k  = _mm256_set1_pd(gain * fs);
    const __m256d lo = _mm256_set1_pd(-fs), hi = _mm256_set1_pd(fs - 1.0);
    const __m256d ch = _mm256_set1_pd(fs - 0.5), cl = _mm256_set1_pd(-fs - 0.5);
    size_t clipped = 0, i = 0;

    for (; i + 4 <= n; i += 4) {
        const __m256d v = _mm256_mul_pd(_mm256_loadu_pd(x + i), k);
        const __m256d c = _mm256_or_pd(_mm256_cmp_pd(v, ch, _CMP_GE_OQ),
                                       _mm256_cmp_pd(v, cl, _CMP_LT_OQ));
        clipped += (size_t)_mm_popcnt_u32((unsigned)_mm256_movemask_pd(c));
        _mm_storeu_si128((__m128i *)(y + i),
                         _mm256_cvtpd_epi32(_mm256_min_pd(_mm256_max_pd(v, lo), hi)));
    }
    return clipped + fu_pcm_to_int_ref(x + i, y + i, n - i, gain, bits);
}

double fu_pcm_peak_ref(const double *x, size_t n)
{
    double m = 0.0;
    for (size_t i = 0; i < n; i++)
        m = fmax(m, fabs(x[i]));
    return m;
}

double fu_pcm_peak(const double *x, size_t n)
{
    const __m256d absmask = _mm256_castsi256_pd(_mm256_set1_epi64x(INT64_MAX));
    __m256d m0 = _mm256_setzero_pd(), m1 = m0, m2 = m0, m3 = m0;
    size_t i = 0;

    for (; i + 16 <= n; i += 16) {
        m0 = _mm256_max_pd(m0, _mm256_and_pd(absmask, _mm256_loadu_pd(x + i)));
        m1 = _mm256_max_pd(m1, _mm256_and_pd(absmask, _mm256_loadu_pd(x + i + 4)));
        m2 = _mm256_max_pd(m2, _mm256_and_pd(absmask, _mm256_loadu_pd(x + i + 8)));
        m3 = _mm256_max_pd(m3, _mm256_and_pd(absmask, _mm256_loadu_pd(x + i + 12)));
    }
    m0 = _mm256_max_pd(_mm256_max_pd(m0, m1), _mm256_max_pd(m2, m3));
    __m128d m = _mm_max_pd(_mm256_castpd256_pd128(m0), _mm256_extractf128_pd(m0, 1));
    m = _mm_max_sd(m, _mm_unpackhi_pd(m, m));
    return fmax(_mm_cvtsd_f64(m), fu_pcm_peak_ref(x + i, n - i));
}
