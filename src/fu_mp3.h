/*
 * fu_mp3.h -- MPEG audio (MP3) reader on libmpg123.
 *
 * Purpose : decode an MP3 to planar doubles straight from the decoder's
 *           floating-point output: no 16-bit step, and peaks above full scale
 *           (which lossy decoding produces) are kept, not clipped.
 * Contract: same shape as the FLAC reader. info.bps = 0 (a lossy stream has no
 *           bit depth). info.frames is exact after a full scan of the file and
 *           already excludes encoder delay and padding when the file carries
 *           that information (gapless). ID3v2 text and pictures (or ID3v1 text)
 *           are exposed as fu_tags.
 */
#ifndef FU_MP3_H
#define FU_MP3_H

#include "fu_audio.h"

#include <stddef.h>

typedef struct fu_mp3_in fu_mp3_in;

fu_mp3_in     *fu_mp3_open(const char *path, fu_format *info, char *err, size_t errlen);
long           fu_mp3_read(fu_mp3_in *in, const double **planes);
const fu_tags *fu_mp3_tags(const fu_mp3_in *in);
void           fu_mp3_close(fu_mp3_in *in);

#endif
