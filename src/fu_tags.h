/*
 * fu_tags.h -- tag and picture store shared by the readers and the FLAC writer.
 *
 * Contract: blocks are owned by the fu_tags and freed by fu_tags_free.
 *           fu_tags_comment appends NAME=value to the (single) VORBIS_COMMENT
 *           block, creating it on first use; value need not be terminated
 *           (len bytes are used). fu_tags_picture appends a PICTURE block and
 *           copies the image. All functions return 0 on success, -1 when the
 *           entry was refused (illegal name, invalid UTF-8, store full).
 */
#ifndef FU_TAGS_H
#define FU_TAGS_H

#include "fu_audio.h"

#include <FLAC/metadata.h>
#include <stddef.h>

enum { FU_TAGS_MAX = 64 };

struct fu_tags {
    FLAC__StreamMetadata *block[FU_TAGS_MAX];
    size_t                n;
};

int  fu_tags_adopt(fu_tags *t, FLAC__StreamMetadata *block);
int  fu_tags_clone(fu_tags *dst, const fu_tags *src);
int  fu_tags_comment(fu_tags *t, const char *name, const char *value, size_t len);
int  fu_tags_picture(fu_tags *t, int type, const char *mime, const char *description,
                     const unsigned char *data, size_t size);
void fu_tags_free(fu_tags *t);

#endif
