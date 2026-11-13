/*
 * fu_pcm.h -- the boundary between double-precision signal and integer PCM.
 *
 * Purpose : scale + round + saturate doubles to `bits`-wide samples (stored in
 *           int32), and find the peak magnitude of a block.
 * Contract: full scale is +-1.0 = +-2^(bits-1), bits in 2..32. Rounding is to
 *           nearest, ties to even. A sample counts as clipped when its rounded
 *           value falls outside [-2^(bits-1), 2^(bits-1) - 1]; it is then
 *           saturated. The *_ref functions are the plain scalar forms; the
 *           others are the AVX2 forms and return identical results (checked by
 *           pcmcheck).
 */
#ifndef FU_PCM_H
#define FU_PCM_H

#include <stddef.h>
#include <stdint.h>

size_t fu_pcm_to_int(const double *x, int32_t *y, size_t n, double gain, unsigned bits);
size_t fu_pcm_to_int_ref(const double *x, int32_t *y, size_t n, double gain, unsigned bits);

double fu_pcm_peak(const double *x, size_t n);
double fu_pcm_peak_ref(const double *x, size_t n);

/* Largest magnitude that still rounds inside the range of `bits`. */
double fu_pcm_limit(unsigned bits);

#endif
