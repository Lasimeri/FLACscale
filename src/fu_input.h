/*
 * fu_input.h -- one reader for every supported source file.
 *
 * Purpose : open a file by what it contains, not by its name:
 *             "fLaC" (after an optional ID3v2 tag)   -> FLAC   (libFLAC)
 *             "RIFF" / "RF64" / "FORM" / "caff"      -> WAV, AIFF, CAF (libsndfile)
 *             anything else                          -> MP3    (libmpg123)
 *           Everything downstream sees planar doubles at +-1.0.
 * Contract: fu_input_read returns frames, 0 at the end, -1 on a decode error;
 *           planes stay valid until the next call. fu_input_kind returns the
 *           container name ("FLAC", "WAV", "AIFF", "MP3", ...).
 */
#ifndef FU_INPUT_H
#define FU_INPUT_H

#include "fu_audio.h"

#include <stddef.h>

typedef struct fu_input fu_input;

fu_input      *fu_input_open(const char *path, fu_format *info, char *err, size_t errlen);
long           fu_input_read(fu_input *in, const double **planes);
const fu_tags *fu_input_tags(const fu_input *in);
const char    *fu_input_kind(const fu_input *in);
void           fu_input_close(fu_input *in);

#endif
