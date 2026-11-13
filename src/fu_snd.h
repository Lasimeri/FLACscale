/*
 * fu_snd.h -- uncompressed-container reader (WAV, RF64, AIFF, CAF) on libsndfile.
 *
 * Contract: same shape as the FLAC and MP3 readers: planar doubles at +-1.0
 *           (integer PCM is divided by 2^(bps-1); float files pass through),
 *           planes valid until the next call, frames / 0 at end / -1 on error.
 *           info.bps is 8, 16, 24 or 32 for integer PCM and 0 for float data.
 *           Text chunks (title, artist, album, date, track, genre) become tags.
 *           fu_snd_kind names the container ("WAV", "AIFF", ...).
 */
#ifndef FU_SND_H
#define FU_SND_H

#include "fu_audio.h"

#include <stddef.h>

typedef struct fu_snd_in fu_snd_in;

fu_snd_in     *fu_snd_open(const char *path, fu_format *info, char *err, size_t errlen);
long           fu_snd_read(fu_snd_in *in, const double **planes);
const fu_tags *fu_snd_tags(const fu_snd_in *in);
const char    *fu_snd_kind(const fu_snd_in *in);
void           fu_snd_close(fu_snd_in *in);

#endif
