/*
 * fu_resample.h -- rational polyphase resampler, the carrier of FLACscale.
 *
 * Purpose : move a planar double-precision signal from rate_in to rate_out
 *           (rate_out >= rate_in) with a Kaiser-windowed sinc, exactly time
 *           aligned (output n sits at input time n * M / L, zero delay).
 * Shape   : a chain of at most two stages.
 *             rate_out >= 2 * rate_in : steep 2x stage, then a relaxed stage
 *                                       for the remaining ratio (if any)
 *             rate_out <  2 * rate_in : one steep stage
 * Contract: frequencies handed to fu_design / fu_stage_init are in cycles per
 *           stage INPUT sample (0.5 = the stage's input Nyquist).
 *           fu_*_run consumes all n_in frames; out must hold fu_*_bound(n_in).
 *           fu_*_flush is called once, after the last run.
 *           Total output frames for n input frames = ceil(n * L / M).
 */
#ifndef FU_RESAMPLE_H
#define FU_RESAMPLE_H

#include "fu_audio.h"

#include <stddef.h>
#include <stdint.h>

typedef struct fu_stage {
    uint32_t L, M, N;          /* up factor, down factor, taps per phase (N % 16 == 0) */
    uint32_t nch;
    uint32_t step;             /* M / L: whole input frames consumed per output */
    uint32_t row;              /* visit-order row of the next output */
    size_t   w;                /* window start of the next output inside buf */
    size_t   fill;             /* valid frames in buf */
    size_t   cap;              /* allocated frames in buf */
    double  *tab;              /* L rows of N taps, stored in the order they are visited */
    uint8_t *adv;              /* per row: extra input frame consumed after this row (0 or 1) */
    double  *buf[FU_MAX_CH];   /* per channel: carried history + new input */
} fu_stage;

typedef struct fu_quality {
    uint32_t taps_steep;       /* taps per phase of the steep stage (default 256) */
    uint32_t taps_relaxed;     /* taps per phase of the relaxed stage (default 32) */
    double   pass;             /* passband edge as a fraction of the source Nyquist (default 0.907) */
} fu_quality;

typedef struct fu_chain {
    uint32_t nstage, nch;
    uint32_t rate_in, rate_out;
    fu_stage st[2];
    double  *mid[FU_MAX_CH];   /* output of stage 0 when two stages run */
    size_t   mid_cap;
} fu_chain;

/* Kaiser attenuation (dB) the design formula assigns to N taps and this transition. */
double fu_design_atten(uint32_t N, double fpass, double fstop);

/* Prototype low-pass at L times the input rate: h[0 .. L*N-1], DC gain L,
 * exactly symmetric about index L*N/2 (h[0] = 0). */
void fu_design(uint32_t L, uint32_t N, double fpass, double fstop, double *h);

int    fu_stage_init(fu_stage *s, uint32_t L, uint32_t M, uint32_t N,
                     double fpass, double fstop, uint32_t nch);
void   fu_stage_free(fu_stage *s);
size_t fu_stage_bound(const fu_stage *s, size_t n_in);
size_t fu_stage_run(fu_stage *s, const double *const in[], size_t n_in,
                    double *const out[], size_t out_cap);
size_t fu_stage_flush(fu_stage *s, double *const out[], size_t out_cap);

fu_quality fu_quality_default(void);
int    fu_chain_init(fu_chain *c, uint32_t rate_in, uint32_t rate_out, uint32_t nch,
                     const fu_quality *q);
void   fu_chain_free(fu_chain *c);
size_t fu_chain_bound(const fu_chain *c, size_t n_in);
size_t fu_chain_run(fu_chain *c, const double *const in[], size_t n_in,
                    double *const out[], size_t out_cap);
size_t fu_chain_flush(fu_chain *c, double *const out[], size_t out_cap);

#endif
