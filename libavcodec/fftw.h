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

#ifndef AVCODEC_FFTW_H
#define AVCODEC_FFTW_H

/**
 * @file
 * @ingroup libavc
 * FFTW based FFT functions header
 */

#include <fftw3.h>

typedef fftwf_complex FFTWComplex;
typedef fftwf_plan FFTWPlan;

typedef struct FFTWContext {
    /**
     * Length of the FFT. Note that the length need not be a power of 2; powers
     * of 2 are most efficient though.
     */
    int n;
    /**
     * if 0 perform forward transform, if 1 perform reverse transform.
     */
    int inverse;
    /**
     * An opaque internal pointer containing plan information, essentially target
     * optimization information for the specified FFT parameters such as n.
     * NOTE: If one calls fft_calc on the exact same array repeatedly, or on the
     * same length with proper alignment, and with the same flags, there is no
     * issue. However, if these conditions are not met, it is necessary to
     * recompute plans. Note that this is often cheap. For precise information on
     * this, see the fftw docs: http://www.fftw.org/fftw3.pdf, or
     * http://www.fftw.org/faq/section3.html#planperarray.
     */
    FFTWPlan plan;
    /**
     * Do a complex, in-place FFT with the parameters defined in ff_fftw_init().
     * No 1.0/sqrt(n) normalization is done.
     */
    void (*fft_calc)(struct FFTWContext *s, FFTWComplex *z);
    /* TODO: add dct, rdft, etc */
} FFTWContext;

/**
 * Set up a complex, in-place FFT.
 * @param s         pointer to an FFTWContext
 * @param n         length of the input array
 * @param inverse   if 0 perform the forward transform, if 1 perform the inverse
 * @return          0 on success, negative AVERROR code otherwise
 */
int ff_fftw_init(FFTWContext *s, int n, int inverse);

/**
 * Clean up an FFTWContext.
 * @param s pointer to an FFTWContext
 */
void ff_fftw_deinit(FFTWContext *s);


#endif /* AVCODEC_FFTW_H */
