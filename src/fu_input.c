/*
 * fu_input.c -- content-sniffing reader front (see fu_input.h).
 */
#include "fu_input.h"
#include "fu_flac.h"
#include "fu_mp3.h"
#include "fu_snd.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum { K_UNREADABLE = -1, K_MP3, K_FLAC, K_SND } kind;

struct fu_input {
    fu_flac_in *flac;
    fu_mp3_in  *mp3;
    fu_snd_in  *snd;
};

static kind sniff(const char *path)
{
    unsigned char h[10];
    FILE *f = fopen(path, "rb");
    if (!f)
        return K_UNREADABLE;
    kind k = K_MP3;
    if (fread(h, 1, 10, f) == 10) {
        if (!memcmp(h, "RIFF", 4) || !memcmp(h, "RF64", 4) || !memcmp(h, "FORM", 4) ||
            !memcmp(h, "caff", 4)) {
            k = K_SND;
        } else {
            long skip = 0;
            if (!memcmp(h, "ID3", 3))   /* ID3v2: 10-byte header + 28-bit syncsafe size */
                skip = 10 + (((long)h[6] & 127) << 21 | ((long)h[7] & 127) << 14 |
                             ((long)h[8] & 127) << 7 | ((long)h[9] & 127));
            if (fseek(f, skip, SEEK_SET) == 0 && fread(h, 1, 4, f) == 4 && !memcmp(h, "fLaC", 4))
                k = K_FLAC;
        }
    }
    fclose(f);
    return k;
}

fu_input *fu_input_open(const char *path, fu_format *info, char *err, size_t errlen)
{
    const kind k = sniff(path);
    fu_input *in = calloc(1, sizeof *in);
    if (k == K_UNREADABLE || !in) {
        snprintf(err, errlen, "%s: %s", path, in ? "cannot read the file; check the path"
                                                 : "out of memory");
        free(in);
        return NULL;
    }
    memset(info, 0, sizeof *info);
    if (k == K_FLAC)
        in->flac = fu_flac_open(path, info, err, errlen);
    else if (k == K_SND)
        in->snd = fu_snd_open(path, info, err, errlen);
    else
        in->mp3 = fu_mp3_open(path, info, err, errlen);
    if (!in->flac && !in->snd && !in->mp3) {
        free(in);
        return NULL;
    }
    return in;
}

long fu_input_read(fu_input *in, const double **planes)
{
    return in->flac ? fu_flac_read(in->flac, planes)
         : in->snd  ? fu_snd_read(in->snd, planes)
                    : fu_mp3_read(in->mp3, planes);
}

const fu_tags *fu_input_tags(const fu_input *in)
{
    return in->flac ? fu_flac_tags(in->flac)
         : in->snd  ? fu_snd_tags(in->snd)
                    : fu_mp3_tags(in->mp3);
}

const char *fu_input_kind(const fu_input *in)
{
    return in->flac ? "FLAC" : in->snd ? fu_snd_kind(in->snd) : "MP3";
}

void fu_input_close(fu_input *in)
{
    if (!in)
        return;
    fu_flac_close(in->flac);
    fu_snd_close(in->snd);
    fu_mp3_close(in->mp3);
    free(in);
}
