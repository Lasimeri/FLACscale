/*
 * fu_tags.c -- tag and picture store (see fu_tags.h).
 */
#include "fu_tags.h"

#include <stdlib.h>
#include <string.h>

int fu_tags_adopt(fu_tags *t, FLAC__StreamMetadata *block)
{
    if (!block || t->n >= FU_TAGS_MAX) {
        if (block)
            FLAC__metadata_object_delete(block);
        return -1;
    }
    t->block[t->n++] = block;
    return 0;
}

int fu_tags_clone(fu_tags *dst, const fu_tags *src)
{
    for (size_t i = 0; src && i < src->n; i++)
        if (fu_tags_adopt(dst, FLAC__metadata_object_clone(src->block[i])))
            return -1;
    return 0;
}

static FLAC__StreamMetadata *comment_block(fu_tags *t)
{
    for (size_t i = 0; i < t->n; i++)
        if (t->block[i]->type == FLAC__METADATA_TYPE_VORBIS_COMMENT)
            return t->block[i];
    FLAC__StreamMetadata *vc = FLAC__metadata_object_new(FLAC__METADATA_TYPE_VORBIS_COMMENT);
    return fu_tags_adopt(t, vc) ? NULL : vc;
}

int fu_tags_comment(fu_tags *t, const char *name, const char *value, size_t len)
{
    FLAC__StreamMetadata *vc = comment_block(t);
    char *text = malloc(len + 1);
    if (!vc || !text) {
        free(text);
        return -1;
    }
    memcpy(text, value, len);
    text[len] = '\0';

    FLAC__StreamMetadata_VorbisComment_Entry e;
    int rc = -1;
    if (FLAC__metadata_object_vorbiscomment_entry_from_name_value_pair(&e, name, text)) {
        if (FLAC__metadata_object_vorbiscomment_append_comment(vc, e, false))
            rc = 0;
        else
            free(e.entry);
    }
    free(text);
    return rc;
}

int fu_tags_picture(fu_tags *t, int type, const char *mime, const char *description,
                    const unsigned char *data, size_t size)
{
    FLAC__StreamMetadata *pic = FLAC__metadata_object_new(FLAC__METADATA_TYPE_PICTURE);
    if (!pic)
        return -1;
    pic->data.picture.type = (FLAC__StreamMetadata_Picture_Type)type;
    if (!FLAC__metadata_object_picture_set_mime_type(pic, (char *)(uintptr_t)mime, true) ||
        !FLAC__metadata_object_picture_set_description(
            pic, (FLAC__byte *)(uintptr_t)description, true) ||
        !FLAC__metadata_object_picture_set_data(pic, (FLAC__byte *)(uintptr_t)data,
                                                (FLAC__uint32)size, true) ||
        !FLAC__metadata_object_picture_is_legal(pic, NULL)) {
        FLAC__metadata_object_delete(pic);
        return -1;
    }
    return fu_tags_adopt(t, pic);
}

void fu_tags_free(fu_tags *t)
{
    for (size_t i = 0; i < t->n; i++)
        FLAC__metadata_object_delete(t->block[i]);
    t->n = 0;
}
