/*
 * fu_cliplevel.c -- clip-level detection (see fu_cliplevel.h).
 */
#include "fu_cliplevel.h"

#include <math.h>
#include <string.h>

enum { HALF = 32768, MIN_COUNT = 32, RATIO = 8, INSIDE = 7 };

void fu_cliplevel_reset(fu_cliplevel *c)
{
    memset(c->hist, 0, sizeof c->hist);
    c->max = -INFINITY;
    c->min = INFINITY;
}

static int bin_of(double x)
{
    const double v = fmin(fmax(x * HALF, -HALF), HALF);
    return (int)lrint(v) + HALF;
}

void fu_cliplevel_feed(fu_cliplevel *c, const double *x, size_t n)
{
    double mx = c->max, mn = c->min;
    for (size_t i = 0; i < n; i++) {
        c->hist[bin_of(x[i])]++;
        mx = fmax(mx, x[i]);
        mn = fmin(mn, x[i]);
    }
    c->max = mx;
    c->min = mn;
}

/* Mean of the INSIDE bins next to `bin`, walking towards zero (dir = -1 or +1). */
static double inside_mean(const fu_cliplevel *c, int bin, int dir)
{
    double sum = 0.0;
    for (int k = 1; k <= INSIDE; k++) {
        const int b = bin + dir * k;
        sum += (b >= 0 && b < FU_CLIPLEVEL_BINS) ? c->hist[b] : 0.0;
    }
    return sum / INSIDE;
}

fu_clipverdict fu_cliplevel_decide(const fu_cliplevel *c)
{
    fu_clipverdict v = { .hi = INFINITY, .lo = -INFINITY };
    if (!(c->max >= c->min))
        return v;
    const int bh = bin_of(c->max), bl = bin_of(c->min);
    v.at_hi = c->hist[bh];
    v.at_lo = c->hist[bl];
    v.inside_hi = inside_mean(c, bh, -1);
    v.inside_lo = inside_mean(c, bl, +1);
    if (v.at_hi >= MIN_COUNT && (double)v.at_hi >= RATIO * fmax(v.inside_hi, 1.0))
        v.hi = ((double)bh - HALF - 0.5) / HALF;
    if (v.at_lo >= MIN_COUNT && (double)v.at_lo >= RATIO * fmax(v.inside_lo, 1.0))
        v.lo = ((double)bl - HALF + 0.5) / HALF;
    return v;
}
