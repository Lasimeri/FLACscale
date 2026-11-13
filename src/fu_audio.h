/*
 * fu_audio.h -- the two types every reader and writer shares.
 *
 * fu_format : stream shape. bps = 0 means floating-point samples with no
 *             integer bit depth (an MP3 decode, a float WAV). lossy = 1 means
 *             the samples come out of a lossy decoder (MP3).
 * fu_tags   : tags and pictures, held as libFLAC metadata blocks
 *             (VORBIS_COMMENT, PICTURE) whatever the source container was.
 *             Opaque here; defined in fu_tags.h for the code that fills it.
 */
#ifndef FU_AUDIO_H
#define FU_AUDIO_H

#include <stdint.h>

#define FU_MAX_CH 8

typedef struct fu_format {
    uint32_t rate, channels, bps, lossy;
    uint64_t frames;            /* 0 when the source does not state it */
} fu_format;

typedef struct fu_tags fu_tags;

#endif
