/*
 * fu_snd.c -- WAV / AIFF / CAF reader on libsndfile 1.2 (see fu_snd.h).
 *
 * Scaling : sf_readf_double with libsndfile's default normalisation returns
 *           integer PCM divided by 2^(bits-1) and float data unchanged, which
 *           is exactly the +-1.0 convention of the other readers.
 */
#include "fu_snd.h"
#include "fu_tags.h"

#include <sndfile.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { CHUNK = 4096 };      /* frames requested per read */

struct fu_snd_in {
    SNDFILE   *sf;
    fu_format  info;
    fu_tags    tags;
    const char *kind;
    double    *inter;                 /* interleaved frames */
    double    *plane[FU_MAX_CH];
};

static uint32_t depth_of(int format)
{
    switch (format & SF_FORMAT_SUBMASK) {
    case SF_FORMAT_PCM_S8:
    case SF_FORMAT_PCM_U8: return 8;
    case SF_FORMAT_PCM_16: return 16;
    case SF_FORMAT_PCM_24: return 24;
    case SF_FORMAT_PCM_32: return 32;
    case SF_FORMAT_FLOAT:
    case SF_FORMAT_DOUBLE: return 0;
    default:               return UINT32_MAX;      /* compressed subtypes: not handled here */
    }
}

static const char *kind_of(int format)
{
    switch (format & SF_FORMAT_TYPEMASK) {
    case SF_FORMAT_WAV:
    case SF_FORMAT_WAVEX:
    case SF_FORMAT_RF64: return "WAV";
    case SF_FORMAT_AIFF: return "AIFF";
    case SF_FORMAT_CAF:  return "CAF";
    default:             return "PCM";
    }
}

static void read_tags(fu_snd_in *in)
{
    static const struct {
        int id;
        const char *name;
    } map[] = { { SF_STR_TITLE, "TITLE" },  { SF_STR_ARTIST, "ARTIST" }, { SF_STR_ALBUM, "ALBUM" },
                { SF_STR_DATE, "DATE" },    { SF_STR_TRACKNUMBER, "TRACKNUMBER" },
                { SF_STR_GENRE, "GENRE" },  { SF_STR_COMMENT, "COMMENT" } };
    for (size_t i = 0; i < sizeof map / sizeof *map; i++) {
        const char *s = sf_get_string(in->sf, map[i].id);
        if (s && *s)
            fu_tags_comment(&in->tags, map[i].name, s, strlen(s));
    }
}

fu_snd_in *fu_snd_open(const char *path, fu_format *info, char *err, size_t errlen)
{
    SF_INFO si = { 0 };
    fu_snd_in *in = calloc(1, sizeof *in);
    if (!in || !(in->sf = sf_open(path, SFM_READ, &si))) {
        snprintf(err, errlen, "%s: cannot open as WAV/AIFF (%s)", path,
                 in ? sf_strerror(NULL) : "out of memory");
        free(in);
        return NULL;
    }
    const uint32_t bps = depth_of(si.format);
    if (bps == UINT32_MAX || si.channels < 1 || si.channels > FU_MAX_CH || si.samplerate < 1) {
        snprintf(err, errlen, "%s: unsupported audio inside the container (format 0x%x, %d "
                              "channels); only PCM and float are read", path, (unsigned)si.format,
                 si.channels);
        fu_snd_close(in);
        return NULL;
    }
    in->kind = kind_of(si.format);
    in->info = (fu_format){ .rate = (uint32_t)si.samplerate, .channels = (uint32_t)si.channels,
                            .bps = bps, .lossy = 0,
                            .frames = si.frames > 0 ? (uint64_t)si.frames : 0 };
    in->inter = malloc((size_t)CHUNK * (size_t)si.channels * sizeof *in->inter);
    for (int ch = 0; ch < si.channels; ch++)
        in->plane[ch] = malloc(CHUNK * sizeof(double));
    if (!in->inter || !in->plane[si.channels - 1]) {
        snprintf(err, errlen, "%s: out of memory", path);
        fu_snd_close(in);
        return NULL;
    }
    read_tags(in);
    *info = in->info;
    return in;
}

long fu_snd_read(fu_snd_in *in, const double **planes)
{
    const uint32_t nch = in->info.channels;
    const sf_count_t n = sf_readf_double(in->sf, in->inter, CHUNK);
    if (n < 0 || (n == 0 && sf_error(in->sf) != SF_ERR_NO_ERROR))
        return -1;
    for (uint32_t ch = 0; ch < nch; ch++) {
        const double *src = in->inter + ch;
        double *dst = in->plane[ch];
        for (sf_count_t i = 0; i < n; i++)
            dst[i] = src[(size_t)i * nch];
        planes[ch] = dst;
    }
    return (long)n;
}

const fu_tags *fu_snd_tags(const fu_snd_in *in)
{
    return &in->tags;
}

const char *fu_snd_kind(const fu_snd_in *in)
{
    return in->kind;
}

void fu_snd_close(fu_snd_in *in)
{
    if (!in)
        return;
    if (in->sf)
        sf_close(in->sf);
    fu_tags_free(&in->tags);
    free(in->inter);
    for (int ch = 0; ch < FU_MAX_CH; ch++)
        free(in->plane[ch]);
    free(in);
}
