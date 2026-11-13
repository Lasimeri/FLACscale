/*
 * fu_lpc.h -- the linear-prediction core of the Opus concealer, in double.
 *
 * Purpose : port of the analysis that celt_decode_lost() runs before it
 *           conceals a frame (libopus 1.6.1, celt/celt_decoder.c:856-877,
 *           celt/celt_lpc.c:37-143 and :281-374, celt/modes.c:392):
 *             taper both ends with the CELT window -> autocorrelation ->
 *             noise floor on lag 0 -> Gaussian lag window -> Levinson-Durbin.
 * Contract: the predictor is e[n] = x[n] + sum_{k=1..order} lpc[k-1] * x[n-k]
 *           (the celt_fir convention). fu_lpc_cfg carries the three numeric
 *           guards of the Opus code; fu_lpc_cfg_opus() returns Opus's own
 *           values (noise floor 1e-4 = -40 dB, lag window 0.008, bail 0.001 =
 *           stop once the prediction gain reaches 30 dB).
 *           fu_lpc_analyze returns the prediction gain ac[0] / error (>= 1).
 * License : the ported routines derive from libopus (BSD-3-Clause), Copyright
 *           (c) 2009-2010 Xiph.Org Foundation, written by Jean-Marc Valin.
 *           The full notice, per routine, is in NOTICE.md.
 */
#ifndef FU_LPC_H
#define FU_LPC_H

#define FU_LPC_ORDER   24     /* CELT_LPC_ORDER */
#define FU_LPC_WINDOW  1024   /* MAX_PERIOD: samples the analysis looks at */
#define FU_LPC_OVERLAP 120    /* CELT overlap at 48 kHz: taper length per side */

typedef struct fu_lpc_cfg {
    double noise_floor;       /* ac[0] *= 1 + noise_floor */
    double lag_window;        /* ac[i] -= ac[i] * (lag_window * i)^2 */
    double bail;              /* stop the recursion when error <= bail * ac[0]; 0 = never */
} fu_lpc_cfg;

fu_lpc_cfg fu_lpc_cfg_opus(void);

/* ac[0..lag] of x[0..n) after tapering `overlap` samples at each end. */
void fu_autocorr(const double *x, int n, int overlap, int lag, double *ac);

/* Levinson-Durbin on ac[0..order]; lpc[0..order). Returns the final error. */
double fu_levinson(const double *ac, int order, double bail, double *lpc);

/* The whole analysis on x[0..n). n < 2 * order + 2 or a silent block gives lpc = 0. */
double fu_lpc_analyze(const double *x, int n, const fu_lpc_cfg *cfg, double lpc[FU_LPC_ORDER]);

#endif
