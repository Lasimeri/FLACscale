/*
 * declipcheck.c -- gate S3: score the declipper against a known truth.
 *
 * Purpose : take real, unclipped audio; clip it artificially so that a chosen
 *           fraction of the samples is pinned; repair it with fu_declip; and
 *           measure, on the pinned samples only, how close the repair comes to
 *           the original:
 *             SDR = 10*log10( sum(truth^2) / sum((truth - estimate)^2) )
 *           "in" scores the clipped signal itself (estimate = ceiling), "out"
 *           the repaired one. Several settings of the Opus numeric guards and
 *           of the pass count are compared on identical input.
 * Input   : any file fu_input reads. It should not be clipped already
 *           (check with clipscan), or the "truth" is not the truth.
 * Usage   : declipcheck FILE
 */
#include "fu_declip.h"
#include "fu_input.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct variant {
    const char   *name;
    fu_declip_cfg cfg;
} variant;

static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return 1e3 * (double)ts.tv_sec + 1e-6 * (double)ts.tv_nsec;
}

/* Smallest magnitude t (on a 2^-15 grid) with at most `frac` of the samples at or above it. */
static double threshold_for(double *const *ch, uint32_t nch, size_t n, double frac)
{
    enum { BINS = 32769 };
    size_t *hist = calloc(BINS, sizeof *hist);
    if (!hist)
        return 1.0;
    for (uint32_t c = 0; c < nch; c++)
        for (size_t i = 0; i < n; i++) {
            const double m = fabs(ch[c][i]) * 32768.0;
            hist[m >= BINS - 1 ? BINS - 1 : (size_t)m]++;
        }
    const size_t want = (size_t)(frac * (double)n * nch);
    size_t acc = 0;
    int b = BINS - 1;
    while (b > 1 && acc + hist[b] <= want)
        acc += hist[b--];
    free(hist);
    return (double)(b + 1) / 32768.0;
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "usage: declipcheck FILE\n");
        return 2;
    }
    char err[512];
    fu_format fmt;
    fu_input *in = fu_input_open(argv[1], &fmt, err, sizeof err);
    if (!in) {
        fprintf(stderr, "declipcheck: %s\n", err);
        return 1;
    }
    double *ch[FU_MAX_CH] = { 0 };
    size_t n = 0, cap = 0;
    for (;;) {
        const double *planes[FU_MAX_CH];
        const long got = fu_input_read(in, planes);
        if (got < 0) {
            fprintf(stderr, "declipcheck: %s: decode error\n", argv[1]);
            return 1;
        }
        if (got == 0)
            break;
        if (n + (size_t)got > cap) {
            cap = 2 * (n + (size_t)got);
            for (uint32_t c = 0; c < fmt.channels; c++)
                if (!(ch[c] = realloc(ch[c], cap * sizeof(double)))) {
                    fprintf(stderr, "declipcheck: out of memory\n");
                    return 1;
                }
        }
        for (uint32_t c = 0; c < fmt.channels; c++)
            memcpy(ch[c] + n, planes[c], (size_t)got * sizeof(double));
        n += (size_t)got;
    }
    fu_input_close(in);
    double *y = malloc(n * sizeof *y);
    if (!y) {
        fprintf(stderr, "declipcheck: out of memory\n");
        return 1;
    }

    /* Every setting is spelled out; none inherits from the library default. */
    const fu_declip_cfg opus = { .lpc = fu_lpc_cfg_opus(), .passes = 1, .guard = 4.0 };
    variant v[5];
    v[0] = (variant){ "Opus guards (floor -40 dB), 1 pass", opus };
    v[1] = (variant){ "Opus guards (floor -40 dB), 2 passes", opus };
    v[1].cfg.passes = 2;
    v[2] = (variant){ "floor -80 dB, 1 pass  [library default]", opus };
    v[2].cfg.lpc.noise_floor = 1e-8;
    v[3] = (variant){ "floor -80 dB, 2 passes", opus };
    v[3].cfg.lpc.noise_floor = 1e-8;
    v[3].cfg.passes = 2;
    v[4] = (variant){ "floor -80 dB, lag window 0.002, 1 pass", opus };
    v[4].cfg.lpc.noise_floor = 1e-8;
    v[4].cfg.lpc.lag_window = 0.002;

    printf("declipcheck: %s\n  %u Hz, %u ch, %zu frames (%.3f s)\n", argv[1], fmt.rate,
           fmt.channels, n, (double)n / fmt.rate);
    static const double fracs[] = { 0.003, 0.009, 0.03 };
    for (size_t f = 0; f < sizeof fracs / sizeof *fracs; f++) {
        const double t = threshold_for(ch, fmt.channels, n, fracs[f]);
        printf("\n  artificial ceiling %.5f (%.2f dBFS): target %.1f %% of samples pinned\n", t,
               20.0 * log10(t), 100.0 * fracs[f]);
        printf("  %-40s %8s %8s %8s | %9s %8s %7s %8s | %8s\n", "setting", "SDR in", "SDR out",
               "gain dB", "raised", "held", "refused", "peak dB", "ms");
        for (size_t k = 0; k < sizeof v / sizeof *v; k++) {
            fu_declip *d = fu_declip_new(&v[k].cfg);
            if (!d)
                return 1;
            double tt = 0.0, e_in = 0.0, e_out = 0.0, ms = 0.0;
            for (uint32_t c = 0; c < fmt.channels; c++) {
                const double *x = ch[c];
                for (size_t i = 0; i < n; i++)
                    y[i] = x[i] >= t ? t : x[i] <= -t ? -t : x[i];
                const double t0 = now_ms();
                fu_declip_span(d, y, n, 0, n, t, -t, 1);
                ms += now_ms() - t0;
                for (size_t i = 0; i < n; i++) {
                    if (fabs(x[i]) < t)
                        continue;
                    const double clipped = x[i] > 0.0 ? t : -t;
                    tt += x[i] * x[i];
                    e_in += (x[i] - clipped) * (x[i] - clipped);
                    e_out += (x[i] - y[i]) * (x[i] - y[i]);
                }
            }
            const fu_declip_stats *s = fu_declip_stats_get(d);
            const double sdr_in = 10.0 * log10(tt / (e_in + 1e-300));
            const double sdr_out = 10.0 * log10(tt / (e_out + 1e-300));
            printf("  %-40s %8.2f %8.2f %+8.2f | %9llu %8llu %7llu %+8.2f | %8.1f\n", v[k].name,
                   sdr_in, sdr_out, sdr_out - sdr_in, (unsigned long long)s->raised,
                   (unsigned long long)s->held, (unsigned long long)s->refused,
                   20.0 * log10(s->peak + 1e-300), ms);
            fu_declip_free(d);
        }
    }
    free(y);
    for (uint32_t c = 0; c < fmt.channels; c++)
        free(ch[c]);
    return 0;
}
