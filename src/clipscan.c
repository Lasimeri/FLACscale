/*
 * clipscan.c -- measurement tool: how a FLAC file sits against its ceiling.
 *
 * Purpose : per channel, the extreme sample values, how many samples sit at
 *           each of the 8 levels nearest the positive and negative extreme,
 *           and the run-length distribution of samples pinned at the extreme
 *           (a run = consecutive samples within TOL steps of the extreme).
 *           A clipped master shows a spike at the extreme and runs of 2+;
 *           a clean one shows a single isolated maximum.
 * Usage   : clipscan FILE.flac [-t TOL]     (TOL in integer sample steps, default 0)
 */
#include "fu_input.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { NCH = 8, NBUCKET = 10 };

typedef struct chan {
    int32_t *pcm;
    size_t   n, cap;
} chan;

static unsigned bucket(size_t len)
{
    /* 1, 2, 3-4, 5-8, 9-16, 17-32, 33-64, 65-128, 129-256, 257+ */
    unsigned b = 0;
    for (size_t top = 1; len > top && b < NBUCKET - 1; top *= 2)
        b++;
    return b;
}

static void runs(const int32_t *x, size_t n, int32_t level, int32_t tol, int positive,
                 size_t hist[NBUCKET], size_t *samples, size_t *longest)
{
    size_t len = 0;
    for (size_t i = 0; i <= n; i++) {
        const int at = i < n && (positive ? x[i] >= level - tol : x[i] <= level + tol);
        if (at) {
            len++;
        } else if (len) {
            hist[bucket(len)]++;
            *samples += len;
            *longest = len > *longest ? len : *longest;
            len = 0;
        }
    }
}

int main(int argc, char **argv)
{
    const char *path = NULL;
    int32_t tol = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc)
            tol = (int32_t)strtol(argv[++i], NULL, 10);
        else
            path = argv[i];
    }
    if (!path) {
        fprintf(stderr, "usage: clipscan FILE.flac [-t TOL]\n");
        return 2;
    }
    char err[512];
    fu_format info;
    fu_input *in = fu_input_open(path, &info, err, sizeof err);
    if (!in) {
        fprintf(stderr, "clipscan: %s\n", err);
        return 1;
    }
    chan c[NCH] = { 0 };
    const double scale = (double)(1ull << ((info.bps ? info.bps : 16) - 1));   /* float sources: scanned on a 16-bit grid */
    for (;;) {
        const double *planes[NCH];
        const long n = fu_input_read(in, planes);
        if (n < 0) {
            fprintf(stderr, "clipscan: %s: decode error\n", path);
            return 1;
        }
        if (n == 0)
            break;
        for (uint32_t ch = 0; ch < info.channels; ch++) {
            if (c[ch].n + (size_t)n > c[ch].cap) {
                c[ch].cap = 2 * (c[ch].n + (size_t)n);
                c[ch].pcm = realloc(c[ch].pcm, c[ch].cap * sizeof(int32_t));
                if (!c[ch].pcm) {
                    fprintf(stderr, "clipscan: out of memory\n");
                    return 1;
                }
            }
            for (long i = 0; i < n; i++)
                c[ch].pcm[c[ch].n + (size_t)i] = (int32_t)lrint(planes[ch][i] * scale);
            c[ch].n += (size_t)n;
        }
    }
    fu_input_close(in);

    printf("clipscan: %s\n  %u Hz, %u bit, %u ch, %zu frames, run tolerance %d step(s)\n", path,
           info.rate, info.bps, info.channels, c[0].n, tol);
    for (uint32_t ch = 0; ch < info.channels; ch++) {
        const int32_t *x = c[ch].pcm;
        const size_t n = c[ch].n;
        int32_t hi = x[0], lo = x[0];
        for (size_t i = 1; i < n; i++) {
            hi = x[i] > hi ? x[i] : hi;
            lo = x[i] < lo ? x[i] : lo;
        }
        size_t near_hi[8] = { 0 }, near_lo[8] = { 0 };
        for (size_t i = 0; i < n; i++) {
            if (hi - x[i] < 8)
                near_hi[hi - x[i]]++;
            if (x[i] - lo < 8)
                near_lo[x[i] - lo]++;
        }
        printf("\n  channel %u: max %d (%.3f dBFS), min %d (%.3f dBFS)\n", ch, hi,
               20.0 * log10((double)hi / scale), lo, 20.0 * log10(-(double)lo / scale));
        printf("    samples at max-0..max-7:");
        for (int k = 0; k < 8; k++)
            printf(" %zu", near_hi[k]);
        printf("\n    samples at min+0..min+7:");
        for (int k = 0; k < 8; k++)
            printf(" %zu", near_lo[k]);
        printf("\n    %-9s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s | %9s %8s\n", "run length", "1", "2",
               "3-4", "5-8", "9-16", "17-32", "33-64", "65-128", "129-256", "257+", "samples", "longest");
        for (int positive = 1; positive >= 0; positive--) {
            size_t hist[NBUCKET] = { 0 }, samples = 0, longest = 0;
            runs(x, n, positive ? hi : lo, tol, positive, hist, &samples, &longest);
            printf("    %-9s ", positive ? "at max" : "at min");
            for (int b = 0; b < NBUCKET; b++)
                printf(" %8zu", hist[b]);
            printf(" | %9zu %8zu   (%.4f %% of the channel)\n", samples, longest,
                   100.0 * (double)samples / (double)n);
        }
        free(c[ch].pcm);
    }
    return 0;
}
