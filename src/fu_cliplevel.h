/*
 * fu_cliplevel.h -- decide whether a channel is clipped, and at what level.
 *
 * Method  : histogram of the samples on a 16-bit grid plus the exact extreme
 *           values. Hard clipping piles samples onto the extreme level: the
 *           extreme bin then holds many times more samples than the bins just
 *           inside it. A clean recording touches its extreme once or a few
 *           times. Rule: the side is clipped when the extreme bin holds at
 *           least 32 samples and at least 8 times the mean of the 7 bins
 *           inside it.
 * Contract: feed every sample of the channel, then decide. hi / lo return the
 *           inner edge of the piled-up extreme bin (a 24-bit limiter ceiling
 *           wobbles inside one 16-bit step; every sample in that bin is pinned)
 *           to hand to fu_declip, or +-infinity for a
 *           side that is not clipped (so no sample ever tests as pinned).
 */
#ifndef FU_CLIPLEVEL_H
#define FU_CLIPLEVEL_H

#include <stddef.h>
#include <stdint.h>

enum { FU_CLIPLEVEL_BINS = 65537 };

typedef struct fu_cliplevel {
    uint32_t hist[FU_CLIPLEVEL_BINS];
    double   max, min;
} fu_cliplevel;

typedef struct fu_clipverdict {
    double   hi, lo;                /* ceilings, or +-infinity */
    uint64_t at_hi, at_lo;          /* samples in the extreme bins */
    double   inside_hi, inside_lo;  /* mean of the 7 bins inside each extreme */
} fu_clipverdict;

void           fu_cliplevel_reset(fu_cliplevel *c);
void           fu_cliplevel_feed(fu_cliplevel *c, const double *x, size_t n);
fu_clipverdict fu_cliplevel_decide(const fu_cliplevel *c);

#endif
