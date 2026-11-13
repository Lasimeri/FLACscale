/*
 * fu_mp3.c -- MP3 reader on libmpg123 1.33 (see fu_mp3.h).
 *
 * Sources : mpg123.h (portable 64-bit API: mpg123_open64, mpg123_length64).
 * Format  : floating-point output is requested with MPG123_FORCE_FLOAT BEFORE
 *           the file is opened; asking afterwards is accepted and ignored (the
 *           decoder has fixed its format by then and keeps delivering 16-bit
 *           integers). The encoding the decoder reports after opening is then
 *           checked, and only 32- or 64-bit float is accepted.
 * Tags    : ID3v2 text frames are mapped to Vorbis comment names by frame id;
 *           APIC frames become FLAC PICTURE blocks; ID3v1 is used only when
 *           the file has no ID3v2.
 */
#include "fu_mp3.h"
#include "fu_tags.h"

#include <mpg123.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { CHUNK = 4096 };      /* frames requested per read */

struct fu_mp3_in {
    mpg123_handle *mh;
    fu_format      info;
    fu_tags        tags;
    size_t         width;                 /* bytes per decoded sample: 4 (float) or 8 (double) */
    void          *inter;                 /* interleaved decoder output */
    double        *plane[FU_MAX_CH];
};

static const struct {
    char id[5];
    const char *name;
} id3_map[] = {
    { "TIT2", "TITLE" },  { "TPE1", "ARTIST" },      { "TALB", "ALBUM" },
    { "TPE2", "ALBUMARTIST" }, { "TRCK", "TRACKNUMBER" }, { "TPOS", "DISCNUMBER" },
    { "TCON", "GENRE" },  { "TDRC", "DATE" },        { "TYER", "DATE" },
    { "TCOM", "COMPOSER" }, { "TSRC", "ISRC" },      { "TBPM", "BPM" },
};

/* Length of the first string in an mpg123_string (fill counts the terminator). */
static size_t first_len(const mpg123_string *s)
{
    return s && s->p && s->fill ? strnlen(s->p, s->fill) : 0;
}

static size_t trimmed(const char *s, size_t max)
{
    size_t n = strnlen(s, max);
    while (n && s[n - 1] == ' ')
        n--;
    return n;
}

static void read_tags(fu_mp3_in *in)
{
    mpg123_id3v1 *v1 = NULL;
    mpg123_id3v2 *v2 = NULL;
    if (mpg123_id3(in->mh, &v1, &v2) != MPG123_OK)
        return;
    if (v2) {
        for (size_t i = 0; i < v2->texts; i++) {
            const mpg123_text *t = &v2->text[i];
            for (size_t m = 0; m < sizeof id3_map / sizeof *id3_map; m++)
                if (!memcmp(t->id, id3_map[m].id, 4) && first_len(&t->text))
                    fu_tags_comment(&in->tags, id3_map[m].name, t->text.p, first_len(&t->text));
        }
        for (size_t i = 0; i < v2->pictures; i++) {
            const mpg123_picture *p = &v2->picture[i];
            fu_tags_picture(&in->tags, p->type, p->mime_type.p ? p->mime_type.p : "",
                            p->description.p ? p->description.p : "", p->data, p->size);
        }
    } else if (v1) {
        const struct {
            const char *name, *text;
            size_t max;
        } f[] = { { "TITLE", v1->title, 30 }, { "ARTIST", v1->artist, 30 },
                  { "ALBUM", v1->album, 30 }, { "DATE", v1->year, 4 } };
        for (size_t i = 0; i < sizeof f / sizeof *f; i++)
            if (trimmed(f[i].text, f[i].max))
                fu_tags_comment(&in->tags, f[i].name, f[i].text, trimmed(f[i].text, f[i].max));
    }
}

fu_mp3_in *fu_mp3_open(const char *path, fu_format *info, char *err, size_t errlen)
{
    int e = MPG123_OK;
    fu_mp3_in *in = calloc(1, sizeof *in);
    if (!in || !(in->mh = mpg123_new(NULL, &e))) {
        snprintf(err, errlen, "%s: cannot start the MP3 decoder (%s)", path,
                 in ? mpg123_plain_strerror(e) : "out of memory");
        free(in);
        return NULL;
    }
    mpg123_param2(in->mh, MPG123_ADD_FLAGS, MPG123_QUIET | MPG123_PICTURE | MPG123_FORCE_FLOAT, 0.0);

    long rate = 0;
    int channels = 0, encoding = 0;
    if (mpg123_open64(in->mh, path) != MPG123_OK ||
        mpg123_getformat(in->mh, &rate, &channels, &encoding) != MPG123_OK ||
        channels < 1 || channels > 2) {
        snprintf(err, errlen, "%s: not a FLAC file and not decodable as MP3 (%s)", path,
                 mpg123_strerror(in->mh));
        fu_mp3_close(in);
        return NULL;
    }
    in->width = encoding == MPG123_ENC_FLOAT_32 ? 4 : encoding == MPG123_ENC_FLOAT_64 ? 8 : 0;
    /* Pin rate, channels and encoding so no later frame can change them. */
    if (!in->width || mpg123_format_none(in->mh) != MPG123_OK ||
        mpg123_format(in->mh, rate, channels, encoding) != MPG123_OK ||
        mpg123_scan(in->mh) != MPG123_OK) {
        snprintf(err, errlen, "%s: the MP3 decoder would not deliver floating point (encoding "
                              "0x%x: %s)", path, (unsigned)encoding, mpg123_strerror(in->mh));
        fu_mp3_close(in);
        return NULL;
    }
    const int64_t frames = mpg123_length64(in->mh);
    in->info = (fu_format){ .rate = (uint32_t)rate, .channels = (uint32_t)channels, .bps = 0, .lossy = 1,
                            .frames = frames > 0 ? (uint64_t)frames : 0 };
    in->inter = malloc((size_t)CHUNK * (size_t)channels * in->width);
    for (int ch = 0; ch < channels; ch++)
        in->plane[ch] = malloc(CHUNK * sizeof(double));
    if (!in->inter || !in->plane[channels - 1]) {
        snprintf(err, errlen, "%s: out of memory", path);
        fu_mp3_close(in);
        return NULL;
    }
    read_tags(in);
    *info = in->info;
    return in;
}

long fu_mp3_read(fu_mp3_in *in, const double **planes)
{
    const uint32_t nch = in->info.channels;
    size_t done = 0;

    while (done == 0) {
        const int rc = mpg123_read(in->mh, in->inter, (size_t)CHUNK * nch * in->width, &done);
        if (rc == MPG123_DONE && done == 0)
            return 0;
        if (rc != MPG123_OK && rc != MPG123_DONE)
            return -1;
    }
    const size_t n = done / (nch * in->width);
    for (uint32_t ch = 0; ch < nch; ch++) {
        double *dst = in->plane[ch];
        if (in->width == 4) {
            const float *src = (const float *)in->inter + ch;
            for (size_t i = 0; i < n; i++)
                dst[i] = (double)src[i * nch];
        } else {
            const double *src = (const double *)in->inter + ch;
            for (size_t i = 0; i < n; i++)
                dst[i] = src[i * nch];
        }
        planes[ch] = dst;
    }
    return (long)n;
}

const fu_tags *fu_mp3_tags(const fu_mp3_in *in)
{
    return &in->tags;
}

void fu_mp3_close(fu_mp3_in *in)
{
    if (!in)
        return;
    if (in->mh) {
        mpg123_close(in->mh);
        mpg123_delete(in->mh);
    }
    fu_tags_free(&in->tags);
    free(in->inter);
    for (int ch = 0; ch < FU_MAX_CH; ch++)
        free(in->plane[ch]);
    free(in);
}
