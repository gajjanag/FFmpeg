/*
 * Copyright (c) 2003, 2007-11 Matteo Frigo
 * Copyright (c) 2003, 2007-11 Massachusetts Institute of Technology.
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with FFmpeg; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "libavutil/mem.h"
#include "fftw.h"

static void fft_calc(struct FFTWContext *s, FFTWComplex *z)
{
    fftwf_execute_dft(s->plan, z, z);
}

av_cold int ff_fftw_init(FFTWContext *s, int n, int inverse)
{
    FFTWComplex *tmp;
    s->n = n;
    s->inverse = inverse;
    s->fft_calc = fft_calc;
    // create a temporary array for computing the plan
    // note: there is no real advantage to doing this on the first calc run;
    // plan creation is not guaranteed to not touch the array at hand.
    tmp = av_malloc(sizeof(*tmp) * n);
    if (!tmp)
        goto fail;
    else {
        /* FIXME: plan creation is not threadsafe */
        /* TODO: some time generate "wisdom" files offline for common arches
         and fft sizes, dump them in a folder, w/o compression as the one time
         plan init can get expensive for long FFT's (seconds for len ~ 8192). */
        if (inverse)
            s->plan = fftwf_plan_dft_1d(n, tmp, tmp, FFTW_BACKWARD, FFTW_PATIENT);
        else
            s->plan = fftwf_plan_dft_1d(n, tmp, tmp, FFTW_FORWARD, FFTW_PATIENT);
        /* apparently never happens in default fftw configurations; just being safe */
        if (!s->plan)
            goto fail;
    }
    av_free(tmp);
    return 0;
fail:
    av_freep(&tmp);
    return AVERROR(ENOMEM);
}

av_cold void ff_fftw_deinit(FFTWContext *s)
{
    fftwf_destroy_plan(s->plan);
    // TODO: doing this complete cleanup may not be ideal; it will bring FFTW
    // to a clean slate and thus future init's won't benefit from some of the
    // accumulated plan knowledge/"wisdom".
    fftwf_cleanup();
}
