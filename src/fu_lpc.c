/*
 * fu_lpc.c -- LPC analysis ported from the Opus concealer (see fu_lpc.h).
 *
 * Differences from libopus, all deliberate:
 *   - double precision, signal scaled to +-1.0 instead of +-32768 float
 *   - the three numeric guards are parameters (fu_lpc_cfg) so the declipper's
 *     ground-truth harness can measure what each one costs or buys
 *   - the taper shrinks when the block is shorter than two overlaps
 */
#include "fu_lpc.h"

#include <math.h>
#include <string.h>

fu_lpc_cfg fu_lpc_cfg_opus(void)
{
    return (fu_lpc_cfg){ .noise_floor = 1e-4, .lag_window = 0.008, .bail = 0.001 };
}

/* celt/modes.c:392 -- the CELT (Vorbis power-complementary) window. */
static double celt_window(int i, int overlap)
{
    const double s = sin(0.5 * M_PI * (i + 0.5) / overlap);
    return sin(0.5 * M_PI * s * s);
}

/* celt/celt_lpc.c:281 _celt_autocorr, float path. */
void fu_autocorr(const double *x, int n, int overlap, int lag, double *ac)
{
    double xx[FU_LPC_WINDOW];
    n = n > FU_LPC_WINDOW ? FU_LPC_WINDOW : n;
    overlap = 2 * overlap > n ? n / 2 : overlap;

    memcpy(xx, x, (size_t)n * sizeof *xx);
    for (int i = 0; i < overlap; i++) {
        const double w = celt_window(i, overlap);
        xx[i] *= w;
        xx[n - i - 1] *= w;
    }
    for (int k = 0; k <= lag; k++) {
        double d = 0.0;
        for (int i = k; i < n; i++)
            d += xx[i] * xx[i - k];
        ac[k] = d;
    }
}

/* celt/celt_lpc.c:37 _celt_lpc, float path. */
double fu_levinson(const double *ac, int order, double bail, double *lpc)
{
    double error = ac[0];

    memset(lpc, 0, (size_t)order * sizeof *lpc);
    if (!(ac[0] > 1e-30))
        return 0.0;
    for (int i = 0; i < order; i++) {
        double rr = 0.0;
        for (int j = 0; j < i; j++)
            rr += lpc[j] * ac[i - j];
        rr += ac[i + 1];
        const double r = -rr / error;
        lpc[i] = r;
        for (int j = 0; j < (i + 1) >> 1; j++) {
            const double t1 = lpc[j], t2 = lpc[i - 1 - j];
            lpc[j]         = t1 + r * t2;
            lpc[i - 1 - j] = t2 + r * t1;
        }
        error -= r * r * error;
        if (error <= bail * ac[0])
            break;
    }
    return error;
}

/* celt/celt_decoder.c:856-877. */
double fu_lpc_analyze(const double *x, int n, const fu_lpc_cfg *cfg, double lpc[FU_LPC_ORDER])
{
    double ac[FU_LPC_ORDER + 1];

    if (n < 2 * FU_LPC_ORDER + 2) {
        memset(lpc, 0, FU_LPC_ORDER * sizeof *lpc);
        return 1.0;
    }
    fu_autocorr(x, n, FU_LPC_OVERLAP, FU_LPC_ORDER, ac);
    ac[0] *= 1.0 + cfg->noise_floor;
    for (int i = 1; i <= FU_LPC_ORDER; i++)
        ac[i] -= ac[i] * (cfg->lag_window * cfg->lag_window) * i * i;
    const double error = fu_levinson(ac, FU_LPC_ORDER, cfg->bail, lpc);
    return error > 0.0 ? ac[0] / error : 1.0;
}
