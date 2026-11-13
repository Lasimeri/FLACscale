/*
 * fu_flac.h -- FLAC files in and out (libFLAC 1.5).
 *
 * Reader  : fu_flac_open reads STREAMINFO and keeps the file's VORBIS_COMMENT
 *           and PICTURE blocks (fu_flac_tags). fu_flac_read decodes one FLAC
 *           frame and hands back planar doubles scaled to +-1.0
 *           (sample / 2^(bps-1)); the planes stay valid until the next call.
 *           Returns frames decoded, 0 at end of stream, -1 on a decode error.
 * Writer  : fu_flac_create opens an encoder at any rate and 8..32 bits, copies
 *           `tags` (may be NULL), appends `extra_tags` ("NAME=value" strings)
 *           and a seek table with one point per 10 s. fu_flac_write takes
 *           planar int32. fu_flac_finish closes the file and returns 0 when
 *           libFLAC reports success.
 * Errors  : open/create return NULL and write one line to err that names the
 *           file and the cause.
 */
#ifndef FU_FLAC_H
#define FU_FLAC_H

#include "fu_audio.h"

#include <stddef.h>

typedef struct fu_flac_in  fu_flac_in;
typedef struct fu_flac_out fu_flac_out;

fu_flac_in    *fu_flac_open(const char *path, fu_format *info, char *err, size_t errlen);
long           fu_flac_read(fu_flac_in *in, const double **planes);
const fu_tags *fu_flac_tags(const fu_flac_in *in);
void           fu_flac_close(fu_flac_in *in);

fu_flac_out *fu_flac_create(const char *path, const fu_format *fmt, const fu_tags *tags,
                            const char *const *extra_tags, size_t n_extra, uint32_t threads,
                            char *err, size_t errlen);
int          fu_flac_write(fu_flac_out *out, const int32_t *const *planes, size_t frames);
int          fu_flac_finish(fu_flac_out *out);

#endif
