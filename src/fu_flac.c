/*
 * fu_flac.c -- FLAC reader and writer on libFLAC (see fu_flac.h).
 *
 * Sources : FLAC/stream_decoder.h, FLAC/stream_encoder.h, FLAC/metadata.h
 *           (libFLAC 1.5.0, API version 14).
 */
#include "fu_flac.h"
#include "fu_tags.h"

#include <FLAC/stream_decoder.h>
#include <FLAC/stream_encoder.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct fu_flac_in {
    FLAC__StreamDecoder *dec;
    fu_format            info;
    fu_tags              tags;              /* cloned VORBIS_COMMENT + PICTURE blocks */
    double              *plane[FU_MAX_CH];
    size_t               cap;               /* frames allocated per plane */
    size_t               got;               /* frames delivered by the last write callback */
    int                  failed;
};

struct fu_flac_out {
    FLAC__StreamEncoder *enc;
    fu_tags              tags;              /* what the encoder writes: comments, pictures, seek table */
};

/* ---- reader ------------------------------------------------------------ */

static FLAC__StreamDecoderWriteStatus on_write(const FLAC__StreamDecoder *dec,
                                               const FLAC__Frame *frame,
                                               const FLAC__int32 *const buffer[], void *client)
{
    (void)dec;
    fu_flac_in *in = client;
    const size_t n = frame->header.blocksize;
    const uint32_t nch = frame->header.channels;

    if (nch != in->info.channels) {
        in->failed = 1;
        return FLAC__STREAM_DECODER_WRITE_STATUS_ABORT;
    }
    if (n > in->cap) {
        for (uint32_t ch = 0; ch < nch; ch++) {
            double *p = realloc(in->plane[ch], n * sizeof *p);
            if (!p) {
                in->failed = 1;
                return FLAC__STREAM_DECODER_WRITE_STATUS_ABORT;
            }
            in->plane[ch] = p;
        }
        in->cap = n;
    }
    const double scale = 1.0 / (double)(1ull << (in->info.bps - 1));
    for (uint32_t ch = 0; ch < nch; ch++) {
        const FLAC__int32 *src = buffer[ch];
        double *dst = in->plane[ch];
        for (size_t i = 0; i < n; i++)
            dst[i] = (double)src[i] * scale;
    }
    in->got = n;
    return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}

static void on_metadata(const FLAC__StreamDecoder *dec, const FLAC__StreamMetadata *m, void *client)
{
    (void)dec;
    fu_flac_in *in = client;

    if (m->type == FLAC__METADATA_TYPE_STREAMINFO) {
        in->info.rate     = m->data.stream_info.sample_rate;
        in->info.channels = m->data.stream_info.channels;
        in->info.bps      = m->data.stream_info.bits_per_sample;
        in->info.frames   = m->data.stream_info.total_samples;
    } else {
        fu_tags_adopt(&in->tags, FLAC__metadata_object_clone(m));
    }
}

static void on_error(const FLAC__StreamDecoder *dec, FLAC__StreamDecoderErrorStatus status,
                     void *client)
{
    (void)dec;
    (void)status;
    ((fu_flac_in *)client)->failed = 1;
}

fu_flac_in *fu_flac_open(const char *path, fu_format *info, char *err, size_t errlen)
{
    fu_flac_in *in = calloc(1, sizeof *in);
    if (!in || !(in->dec = FLAC__stream_decoder_new())) {
        snprintf(err, errlen, "%s: out of memory", path);
        free(in);
        return NULL;
    }
    FLAC__stream_decoder_set_metadata_respond(in->dec, FLAC__METADATA_TYPE_VORBIS_COMMENT);
    FLAC__stream_decoder_set_metadata_respond(in->dec, FLAC__METADATA_TYPE_PICTURE);

    const FLAC__StreamDecoderInitStatus st =
        FLAC__stream_decoder_init_file(in->dec, path, on_write, on_metadata, on_error, in);
    if (st != FLAC__STREAM_DECODER_INIT_STATUS_OK) {
        snprintf(err, errlen, "%s: cannot open (%s); check the path", path,
                 FLAC__StreamDecoderInitStatusString[st]);
        fu_flac_close(in);
        return NULL;
    }
    if (!FLAC__stream_decoder_process_until_end_of_metadata(in->dec) || in->failed ||
        in->info.rate == 0 || in->info.channels == 0 || in->info.channels > FU_MAX_CH ||
        in->info.bps < 4 || in->info.bps > 32) {
        snprintf(err, errlen, "%s: no usable FLAC STREAMINFO (%s)", path,
                 FLAC__stream_decoder_get_resolved_state_string(in->dec));
        fu_flac_close(in);
        return NULL;
    }
    *info = in->info;
    return in;
}

long fu_flac_read(fu_flac_in *in, const double **planes)
{
    in->got = 0;
    while (in->got == 0) {
        if (FLAC__stream_decoder_get_state(in->dec) == FLAC__STREAM_DECODER_END_OF_STREAM)
            return 0;
        if (!FLAC__stream_decoder_process_single(in->dec) || in->failed)
            return -1;
    }
    for (uint32_t ch = 0; ch < in->info.channels; ch++)
        planes[ch] = in->plane[ch];
    return (long)in->got;
}

const fu_tags *fu_flac_tags(const fu_flac_in *in)
{
    return &in->tags;
}

void fu_flac_close(fu_flac_in *in)
{
    if (!in)
        return;
    if (in->dec) {
        FLAC__stream_decoder_finish(in->dec);
        FLAC__stream_decoder_delete(in->dec);
    }
    fu_tags_free(&in->tags);
    for (int ch = 0; ch < FU_MAX_CH; ch++)
        free(in->plane[ch]);
    free(in);
}

/* ---- writer ------------------------------------------------------------ */

static void out_free(fu_flac_out *out)
{
    if (out->enc)
        FLAC__stream_encoder_delete(out->enc);
    fu_tags_free(&out->tags);
    free(out);
}

/* Source blocks, then the extra comments, then a seek table template. */
static int build_metadata(fu_flac_out *out, const fu_format *fmt, const fu_tags *tags,
                          const char *const *extra_tags, size_t n_extra)
{
    if (fu_tags_clone(&out->tags, tags))
        return -1;
    for (size_t i = 0; i < n_extra; i++) {
        const char *eq = strchr(extra_tags[i], '=');
        char name[64];
        if (!eq || (size_t)(eq - extra_tags[i]) >= sizeof name)
            return -1;
        memcpy(name, extra_tags[i], (size_t)(eq - extra_tags[i]));
        name[eq - extra_tags[i]] = '\0';
        if (fu_tags_comment(&out->tags, name, eq + 1, strlen(eq + 1)))
            return -1;
    }
    if (fmt->frames) {
        FLAC__StreamMetadata *seek = FLAC__metadata_object_new(FLAC__METADATA_TYPE_SEEKTABLE);
        if (!seek ||
            !FLAC__metadata_object_seektable_template_append_spaced_points_by_samples(
                seek, 10 * fmt->rate, fmt->frames) ||
            !FLAC__metadata_object_seektable_template_sort(seek, true)) {
            if (seek)
                FLAC__metadata_object_delete(seek);
            return -1;
        }
        if (fu_tags_adopt(&out->tags, seek))
            return -1;
    }
    return 0;
}

fu_flac_out *fu_flac_create(const char *path, const fu_format *fmt, const fu_tags *tags,
                            const char *const *extra_tags, size_t n_extra, uint32_t threads,
                            char *err, size_t errlen)
{
    fu_flac_out *out = calloc(1, sizeof *out);
    if (!out || !(out->enc = FLAC__stream_encoder_new()) ||
        build_metadata(out, fmt, tags, extra_tags, n_extra)) {
        snprintf(err, errlen, "%s: cannot prepare the encoder metadata (out of memory or an "
                              "illegal tag)", path);
        if (out)
            out_free(out);
        return NULL;
    }
    FLAC__StreamEncoder *e = out->enc;
    FLAC__stream_encoder_set_channels(e, fmt->channels);
    FLAC__stream_encoder_set_bits_per_sample(e, fmt->bps);
    FLAC__stream_encoder_set_sample_rate(e, fmt->rate);
    FLAC__stream_encoder_set_compression_level(e, 5);
    FLAC__stream_encoder_set_total_samples_estimate(e, fmt->frames);
    FLAC__stream_encoder_set_metadata(e, out->tags.block, (uint32_t)out->tags.n);
    /* libFLAC refuses a count above its own maximum and then encodes on ONE
     * thread (measured with libFLAC 1.5.0: 6.53 s at 96 and at 128 threads, the
     * single-thread time, against 0.63 s at 64). Step down to the largest count
     * it accepts. */
    for (uint32_t t = threads; t > 1; t--)
        if (FLAC__stream_encoder_set_num_threads(e, t) == FLAC__STREAM_ENCODER_SET_NUM_THREADS_OK)
            break;

    const FLAC__StreamEncoderInitStatus st = FLAC__stream_encoder_init_file(e, path, NULL, NULL);
    if (st != FLAC__STREAM_ENCODER_INIT_STATUS_OK) {
        snprintf(err, errlen, "%s: cannot create (%s); check that the folder exists and is writable",
                 path, FLAC__StreamEncoderInitStatusString[st]);
        out_free(out);
        return NULL;
    }
    return out;
}

int fu_flac_write(fu_flac_out *out, const int32_t *const *planes, size_t frames)
{
    return FLAC__stream_encoder_process(out->enc, (const FLAC__int32 *const *)planes,
                                        (uint32_t)frames) ? 0 : -1;
}

int fu_flac_finish(fu_flac_out *out)
{
    const int ok = FLAC__stream_encoder_finish(out->enc);
    out_free(out);
    return ok ? 0 : -1;
}
