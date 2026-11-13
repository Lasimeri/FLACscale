/*
 * fu_dequant.h -- back-calculate precision that rounding to a bit depth removed.
 *
 * Idea    : a quantized sample q does not say "the value was q". It says "the
 *           value was somewhere in [q - half, q + half] steps". Those intervals
 *           are hard facts recovered from the file (Thao and Vetterli, IEEE
 *           Trans. Signal Processing 42, 1994: any estimate is improved by
 *           projecting it into them).
 * Method  : 1. Wiener filter in a time-frequency frame. Rounding noise is flat;
 *              music is not. Where a coefficient's local power is far above
 *              the noise it is kept; where the noise stands alone it is
 *              removed. Gain = S / (S + noise).
 *           2. Projection: every output sample is clamped into its interval.
 *              The truth lies in every interval, so this can only move the
 *              estimate towards it.
 * Rejected: the published sparse formulations (minimise the L1 norm of the
 *           coefficients inside the box, Douglas-Rachford; Zaviska, Rajmic,
 *           Mokry, ICASSP 2021) and a spectrum-weighted variant of them were
 *           built and measured first. On dense music at 8 and 16 bits both
 *           came out WORSE than the rounded signal (-5.4 to -5.9 dB and -2.0
 *           to -5.0 dB): a penalty that prefers "smaller" slides every sample
 *           to the edge of its interval nearest zero. The numbers are in
 *           docs/MEASUREMENTS.md.
 * Contract: x[0..n) holds quantized samples at full scale +-1.0 (multiples of
 *           `step`) and is replaced by the estimate. Every output sample stays
 *           within `half` steps of its input. The signal is solved in
 *           independent overlapping segments (cross-faded, which keeps the
 *           result consistent) on `threads` threads.
 *           Returns 0, or -1 when memory runs out (x is then untouched).
 * Limits  : the gain is small where the music is loud (it covers the noise)
 *           and grows as the passage gets quieter. `half` = 0.5 and `noise` =
 *           1/12 are exact for plain rounding. A dithered master moved each
 *           sample by up to 1.5 steps before rounding (TPDF): the box must
 *           then be 1.5 wide and the noise 1/4.
 */
#ifndef FU_DEQUANT_H
#define FU_DEQUANT_H

#include <stddef.h>

typedef struct fu_dequant_cfg {
    int    passes;     /* 1 = filter once; more = re-estimate the signal power from the result
                          (default 2: measured +0.02 to +0.19 dB over 1; a third adds <= 0.05) */
    double half;       /* half-width of the consistency interval, in steps */
    double noise;      /* rounding-noise power per sample, in steps^2 */
    int    threads;    /* 0 = one per hardware thread */
} fu_dequant_cfg;

fu_dequant_cfg fu_dequant_cfg_default(void);
int            fu_dequant(double *x, size_t n, double step, const fu_dequant_cfg *cfg);

#endif
