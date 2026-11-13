/*
 * flacdiff.c -- measurement tool: how far apart two audio files are.
 *
 * Purpose : null test. Finds the time offset between A and B (cross-correlation
 *           of a 65536-frame excerpt from the middle, lags -4096..4096), the
 *           least-squares gain g that maps B onto A, and the residual A - g*B
 *           over the whole overlap, per channel.
 * Units   : residual dB = 10*log10(sum((A - g*B)^2) / sum(A^2)); -inf means
 *           identical. "worst" is the largest single-sample difference in
 *           units of A's full scale.
 * Usage   : flacdiff A B        (FLAC or MP3, same sample rate and channels)
 */
#include "fu_input.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct track {
    fu_format fmt;
    double   *ch[FU_MAX_CH];
    size_t    n;
} track;

static int load(const char *path, track *t)
{
    char err[512];
    size_t cap = 0;
    fu_input *in = fu_input_open(path, &t->fmt, err, sizeof err);
    if (!in) {
        fprintf(stderr, "flacdiff: %s\n", err);
        return -1;
    }
    for (;;) {
        const double *planes[FU_MAX_CH];
        const long n = fu_input_read(in, planes);
        if (n < 0) {
            fprintf(stderr, "flacdiff: %s: decode error\n", path);
            return -1;
        }
        if (n == 0)
            break;
        if (t->n + (size_t)n > cap) {
            cap = 2 * (t->n + (size_t)n);
            for (uint32_t c = 0; c < t->fmt.channels; c++) {
                t->ch[c] = realloc(t->ch[c], cap * sizeof(double));
                if (!t->ch[c]) {
                    fprintf(stderr, "flacdiff: out of memory\n");
                    return -1;
                }
            }
        }
        for (uint32_t c = 0; c < t->fmt.channels; c++)
            memcpy(t->ch[c] + t->n, planes[c], (size_t)n * sizeof(double));
        t->n += (size_t)n;
    }
    fu_input_close(in);
    return 0;
}

int main(int argc, char **argv)
{
    enum { WIN = 65536, MAXLAG = 4096 };
    track a = { 0 }, b = { 0 };
    if (argc != 3) {
        fprintf(stderr, "usage: flacdiff A B\n");
        return 2;
    }
    if (load(argv[1], &a) || load(argv[2], &b))
        return 1;
    if (a.fmt.rate != b.fmt.rate || a.fmt.channels != b.fmt.channels) {
        fprintf(stderr, "flacdiff: formats differ (%u Hz %u ch vs %u Hz %u ch); resample one first\n",
                a.fmt.rate, a.fmt.channels, b.fmt.rate, b.fmt.channels);
        return 1;
    }
    const size_t nmin = a.n < b.n ? a.n : b.n;
    if (nmin < (size_t)WIN + 2 * MAXLAG) {
        fprintf(stderr, "flacdiff: files are shorter than %d frames\n", WIN + 2 * MAXLAG);
        return 1;
    }
    /* B[i + lag] lines up with A[i]. */
    const size_t mid = nmin / 2 - WIN / 2;
    long lag = 0;
    double best = -1.0;
    for (long l = -MAXLAG; l <= MAXLAG; l++) {
        double s = 0.0;
        for (size_t i = 0; i < WIN; i++)
            s += a.ch[0][mid + i] * b.ch[0][(size_t)((long)(mid + i) + l)];
        if (fabs(s) > best) {
            best = fabs(s);
            lag = l;
        }
    }
    const size_t a0 = lag < 0 ? (size_t)-lag : 0, b0 = lag > 0 ? (size_t)lag : 0;
    const size_t n = (a.n - a0) < (b.n - b0) ? (a.n - a0) : (b.n - b0);

    printf("flacdiff: A = %s\n          B = %s\n", argv[1], argv[2]);
    printf("  A %zu frames, B %zu frames, %u Hz; B leads A by %ld frames; %zu frames compared\n",
           a.n, b.n, a.fmt.rate, lag, n);
    printf("  %-8s %12s %14s %12s\n", "channel", "gain B->A dB", "residual dB", "worst");
    for (uint32_t c = 0; c < a.fmt.channels; c++) {
        const double *x = a.ch[c] + a0, *y = b.ch[c] + b0;
        double xy = 0.0, yy = 0.0, xx = 0.0;
        for (size_t i = 0; i < n; i++) {
            xy += x[i] * y[i];
            yy += y[i] * y[i];
            xx += x[i] * x[i];
        }
        const double g = yy > 0.0 ? xy / yy : 0.0;
        double rr = 0.0, worst = 0.0;
        for (size_t i = 0; i < n; i++) {
            const double r = x[i] - g * y[i];
            rr += r * r;
            worst = fabs(r) > worst ? fabs(r) : worst;
        }
        printf("  %-8u %12.4f %14.2f %12.3e\n", c, 20.0 * log10(fabs(g) + 1e-300),
               10.0 * log10(rr / (xx + 1e-300) + 1e-300), worst);
    }
    return 0;
}
