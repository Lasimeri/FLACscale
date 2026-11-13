/*
 * fu_resample_ref.h -- the plain forms the production kernel is measured against.
 *
 * Purpose : (1) fu_ref_scalar: the resampling sum written straight from the
 *               formula, one channel, no tables, no buffers. It is the oracle
 *               for correctness (values, alignment, output count).
 *           (2) fu_ref_phase2: the same AVX2 stereo dot product as production,
 *               but with the table stored in phase order and indexed by
 *               p = (n*M + c) mod L. It isolates what the visit-order layout
 *               buys.
 * Contract: both produce ceil(n_in * L / M) frames into y; h is the prototype
 *           from fu_design (L*N taps). Used by rsbench only.
 */
#ifndef FU_RESAMPLE_REF_H
#define FU_RESAMPLE_REF_H

#include <stddef.h>
#include <stdint.h>

size_t fu_ref_count(size_t n_in, uint32_t L, uint32_t M);

size_t fu_ref_scalar(uint32_t L, uint32_t M, uint32_t N, const double *h,
                     const double *x, size_t n_in, double *y);

size_t fu_ref_phase2(uint32_t L, uint32_t M, uint32_t N, const double *h,
                     const double *xl, const double *xr, size_t n_in,
                     double *yl, double *yr);

#endif
