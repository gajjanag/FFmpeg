/*
 * AAC encoder utilities
 * Copyright (C) 2015 Rostislav Pehlivanov
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/**
 * @file
 * AAC encoder utilities
 * @author Rostislav Pehlivanov ( atomnuker gmail com )
 */

#ifndef AVCODEC_AACENC_UTILS_H
#define AVCODEC_AACENC_UTILS_H

#include "libavutil/intfloat.h"
#include "aac.h"
#include "aacenctab.h"
#include "aactab.h"

#define ROUND_STANDARD 0.4054f
#define ROUND_TO_ZERO 0.1054f
#define C_QUANT 0.4054f

static inline void abs_pow34_v(float *out, const float *in, const int size)
{
    int i;
    for (i = 0; i < size; i++) {
        float a = fabsf(in[i]);
        out[i] = sqrtf(a * sqrtf(a));
    }
}

static inline float pos_pow34(float a)
{
    return sqrtf(a * sqrtf(a));
}

/**
 * Quantize one coefficient.
 * @return absolute value of the quantized coefficient
 * @see 3GPP TS26.403 5.6.2 "Scalefactor determination"
 */
static inline int quant(float coef, const float Q, const float rounding)
{
    float a = coef * Q;
    return sqrtf(a * sqrtf(a)) + rounding;
}

static inline void quantize_bands(int *out, const float *in, const float *scaled,
                                  int size, float Q34, int is_signed, int maxval,
                                  const float rounding)
{
    int i;
    for (i = 0; i < size; i++) {
        float qc = scaled[i] * Q34;
        out[i] = (int)FFMIN(qc + rounding, (float)maxval);
        if (is_signed && in[i] < 0.0f) {
            out[i] = -out[i];
        }
    }
}

static inline float find_max_val(int group_len, int swb_size, const float *scaled)
{
    float maxval = 0.0f;
    int w2, i;
    for (w2 = 0; w2 < group_len; w2++) {
        for (i = 0; i < swb_size; i++) {
            maxval = FFMAX(maxval, scaled[w2*128+i]);
        }
    }
    return maxval;
}

static inline int find_min_book(float maxval, int sf)
{
    float Q34 = ff_aac_pow34sf_tab[POW_SF2_ZERO - sf + SCALE_ONE_POS - SCALE_DIV_512];
    int qmaxval, cb;
    qmaxval = maxval * Q34 + C_QUANT;
    if (qmaxval >= (FF_ARRAY_ELEMS(aac_maxval_cb)))
        cb = 11;
    else
        cb = aac_maxval_cb[qmaxval];
    return cb;
}

static const float
bp[] = {1.0, 1.5,},
dp_h[] = { 0.0, 5.84960938e-01,}, /* 0x3f15c000 */
dp_l[] = { 0.0, 1.56322085e-06,}, /* 0x35d1cfdc */
zero    =  0.0,
one	=  1.0,
two	=  2.0,
two24	=  16777216.0,	/* 0x4b800000 */
huge	=  1.0e30,
tiny    =  1.0e-30,
	/* poly coefs for (3/2)*(log(x)-2s-2/3*s**3 */
L1  =  6.0000002384e-01, /* 0x3f19999a */
L2  =  4.2857143283e-01, /* 0x3edb6db7 */
L3  =  3.3333334327e-01, /* 0x3eaaaaab */
L4  =  2.7272811532e-01, /* 0x3e8ba305 */
L5  =  2.3066075146e-01, /* 0x3e6c3255 */
L6  =  2.0697501302e-01, /* 0x3e53f142 */
P1   =  1.6666667163e-01, /* 0x3e2aaaab */
P2   = -2.7777778450e-03, /* 0xbb360b61 */
P3   =  6.6137559770e-05, /* 0x388ab355 */
P4   = -1.6533901999e-06, /* 0xb5ddea0e */
P5   =  4.1381369442e-08, /* 0x3331bb4c */
lg2  =  6.9314718246e-01, /* 0x3f317218 */
lg2_h  =  6.93145752e-01, /* 0x3f317200 */
lg2_l  =  1.42860654e-06, /* 0x35bfbe8c */
ovt =  4.2995665694e-08, /* -(128-log2(ovfl+.5ulp)) */
cp    =  9.6179670095e-01, /* 0x3f76384f =2/(3ln2) */
cp_h  =  9.6191406250e-01, /* 0x3f764000 =12b cp */
cp_l  = -1.1736857402e-04, /* 0xb8f623c6 =tail of cp_h */
ivln2    =  1.4426950216e+00, /* 0x3fb8aa3b =1/ln2 */
ivln2_h  =  1.4426879883e+00, /* 0x3fb8aa00 =16b 1/ln2*/
ivln2_l  =  7.0526075433e-06; /* 0x36eca570 =1/ln2 tail*/

#define GET_FLOAT_WORD(word,value)					\
do {                                        \
  union av_intfloat32 data;					\
  data.f = (value);						\
  (word) = data.i;						\
} while (0)

#define SET_FLOAT_WORD(value,word)					\
do {								\
  union av_intfloat32 data;					\
  data.i = (word);						\
  (value) = data.f;						\
} while (0)

static float ff_powf(float x, float y)
{
	float z,ax,z_h,z_l,p_h,p_l;
	float y1,t1,t2,r,s,t,u,v,w;
    float s2,s_h,s_l,t_h,t_l;
	int32_t i,j,k,n;
	int32_t hx,ix,is;

	GET_FLOAT_WORD(hx,x);
	ix = hx;

	ax   = x;

	n = ((uint32_t)hx>>31)-1;

	    n = 0;
	/* take care subnormal number */
	    if(ix<0x00800000)
		{ax *= two24; n -= 24; GET_FLOAT_WORD(ix,ax); }
	    n  += ((ix)>>23)-0x7f;
	    j  = ix&0x007fffff;
	/* determine interval */
	    ix = j|0x3f800000;		/* normalize ix */
	    k=0;/* |x|<sqrt(3/2) */
	    SET_FLOAT_WORD(ax,ix);

	/* compute s = s_h+s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5) */
	    u = ax-bp[k];		/* bp[0]=1.0, bp[1]=1.5 */
	    v = one/(ax+bp[k]);
	    s = u*v;
	    s_h = s;
	    GET_FLOAT_WORD(is,s_h);
	    SET_FLOAT_WORD(s_h,is&0xfffff000);
	/* t_h=ax+bp[k] High */
	    is = ((ix>>1)&0xfffff000)|0x20000000;
	    SET_FLOAT_WORD(t_h,is+0x00400000+(k<<21));
	    t_l = ax - (t_h-bp[k]);
	    s_l = v*((u-s_h*t_h)-s_h*t_l);
	/* compute log(ax) */
	    s2 = s*s;
	    r = s2*s2*(L1+s2*(L2+s2*(L3+s2*(L4+s2*(L5+s2*L6)))));
	    r += s_l*(s_h+s);
	    s2  = s_h*s_h;
	    t_h = (float)3.0+s2+r;
	    GET_FLOAT_WORD(is,t_h);
	    SET_FLOAT_WORD(t_h,is&0xfffff000);
	    t_l = r-((t_h-(float)3.0)-s2);
	/* u+v = s*(1+...) */
	    u = s_h*t_h;
	    v = s_l*t_h+t_l*s;
	/* 2/(3log2)*(s+...) */
	    p_h = u+v;
	    GET_FLOAT_WORD(is,p_h);
	    SET_FLOAT_WORD(p_h,is&0xfffff000);
	    p_l = v-(p_h-u);
	    z_h = cp_h*p_h;		/* cp_h+cp_l = 2/(3*log2) */
	    z_l = cp_l*p_h+p_l*cp+dp_l[k];
	/* log2(ax) = (s+..)*2/(3*log2) = n + dp_h + z_h + z_l */
	    t = (float)n;
	    t1 = (((z_h+z_l)+dp_h[k])+t);
	    GET_FLOAT_WORD(is,t1);
	    SET_FLOAT_WORD(t1,is&0xfffff000);
	    t2 = z_l-(((t1-t)-dp_h[k])-z_h);

    /* split up y into y1+y2 and compute (y1+y2)*(t1+t2) */
	GET_FLOAT_WORD(is,y);
	SET_FLOAT_WORD(y1,is&0xfffff000);
	p_l = (y-y1)*t1+y*t2;
	p_h = y1*t1;
	z = p_l+p_h;
	GET_FLOAT_WORD(j,z);
	if ((j&0x7fffffff)>0x43160000)		/* z <= -150 */
	    return tiny*tiny;			/* underflow */
    else if (j==0xc3160000){			/* z == -150 */
	    if(p_l<=z-p_h) return tiny*tiny;		/* underflow */
	}
    /*
     * compute 2**(p_h+p_l)
     */
	i = j&0x7fffffff;
	k = (i>>23)-0x7f;
	n = 0;
	if(i>0x3f000000) {		/* if |z| > 0.5, set n = [z+0.5] */
	    n = j+(0x00800000>>(k+1));
	    k = ((n&0x7fffffff)>>23)-0x7f;	/* new k for n */
	    SET_FLOAT_WORD(t,n&~(0x007fffff>>k));
	    n = ((n&0x007fffff)|0x00800000)>>(23-k);
	    if(j<0) n = -n;
	    p_h -= t;
	}
	t = p_l+p_h;
	GET_FLOAT_WORD(is,t);
	SET_FLOAT_WORD(t,is&0xffff8000);
	u = t*lg2_h;
	v = (p_l-(t-p_h))*lg2+t*lg2_l;
	z = u+v;
	w = v-(z-u);
	t  = z*z;
	t1  = z - t*(P1+t*(P2+t*(P3+t*(P4+t*P5))));
	r  = (z*t1)/(t1-two)-(w+z*w);
	z  = one-(r-z);
	GET_FLOAT_WORD(j,z);
	j += (n<<23);
	if((j>>23)<=0) z = scalbnf(z,n);	/* subnormal output */
	else SET_FLOAT_WORD(z,j);
	return z;
}

static inline float find_form_factor(int group_len, int swb_size, float thresh,
                                     const float *scaled, float nzslope) {
    const float iswb_size = 1.0f / swb_size;
    const float iswb_sizem1 = 1.0f / (swb_size - 1);
    const float ethresh = thresh;
    float form = 0.0f, weight = 0.0f;
    int w2, i;
    for (w2 = 0; w2 < group_len; w2++) {
        float e = 0.0f, e2 = 0.0f, var = 0.0f, maxval = 0.0f;
        float nzl = 0;
        for (i = 0; i < swb_size; i++) {
            float s = fabsf(scaled[w2*128+i]);
            maxval = FFMAX(maxval, s);
            e += s;
            e2 += s *= s;
            /* We really don't want a hard non-zero-line count, since
             * even below-threshold lines do add up towards band spectral power.
             * So, fall steeply towards zero, but smoothly
             */
            if (s >= ethresh) {
                nzl += 1.0f;
            } else {
                if (nzslope == 2.f)
                    nzl += (s / ethresh) * (s / ethresh);
                else
                    nzl += ff_powf(s / ethresh, nzslope);
                //nzl += expf(logf(s / ethresh) * nzslope);
            }
        }
        if (e2 > thresh) {
            float frm;
            e *= iswb_size;

            /** compute variance */
            for (i = 0; i < swb_size; i++) {
                float d = fabsf(scaled[w2*128+i]) - e;
                var += d*d;
            }
            var = sqrtf(var * iswb_sizem1);

            e2 *= iswb_size;
            frm = e / FFMIN(e+4*var,maxval);
            form += e2 * sqrtf(frm) / FFMAX(0.5f,nzl);
            weight += e2;
        }
    }
    if (weight > 0) {
        return form / weight;
    } else {
        return 1.0f;
    }
}

/** Return the minimum scalefactor where the quantized coef does not clip. */
static inline uint8_t coef2minsf(float coef)
{
    return av_clip_uint8(log2f(coef)*4 - 69 + SCALE_ONE_POS - SCALE_DIV_512);
}

/** Return the maximum scalefactor where the quantized coef is not zero. */
static inline uint8_t coef2maxsf(float coef)
{
    return av_clip_uint8(log2f(coef)*4 +  6 + SCALE_ONE_POS - SCALE_DIV_512);
}

/*
 * Returns the closest possible index to an array of float values, given a value.
 */
static inline int quant_array_idx(const float val, const float *arr, const int num)
{
    int i, index = 0;
    float quant_min_err = INFINITY;
    for (i = 0; i < num; i++) {
        float error = (val - arr[i])*(val - arr[i]);
        if (error < quant_min_err) {
            quant_min_err = error;
            index = i;
        }
    }
    return index;
}

/**
 * approximates exp10f(-3.0f*(0.5f + 0.5f * cosf(FFMIN(b,15.5f) / 15.5f)))
 */
static av_always_inline float bval2bmax(float b)
{
    return 0.001f + 0.0035f * (b*b*b) / (15.5f*15.5f*15.5f);
}

/*
 * Compute a nextband map to be used with SF delta constraint utilities.
 * The nextband array should contain 128 elements, and positions that don't
 * map to valid, nonzero bands of the form w*16+g (with w being the initial
 * window of the window group, only) are left indetermined.
 */
static inline void ff_init_nextband_map(const SingleChannelElement *sce, uint8_t *nextband)
{
    unsigned char prevband = 0;
    int w, g;
    /** Just a safe default */
    for (g = 0; g < 128; g++)
        nextband[g] = g;

    /** Now really navigate the nonzero band chain */
    for (w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]) {
        for (g = 0; g < sce->ics.num_swb; g++) {
            if (!sce->zeroes[w*16+g] && sce->band_type[w*16+g] < RESERVED_BT)
                prevband = nextband[prevband] = w*16+g;
        }
    }
    nextband[prevband] = prevband; /* terminate */
}

/*
 * Updates nextband to reflect a removed band (equivalent to
 * calling ff_init_nextband_map after marking a band as zero)
 */
static inline void ff_nextband_remove(uint8_t *nextband, int prevband, int band)
{
    nextband[prevband] = nextband[band];
}

/*
 * Checks whether the specified band could be removed without inducing
 * scalefactor delta that violates SF delta encoding constraints.
 * prev_sf has to be the scalefactor of the previous nonzero, nonspecial
 * band, in encoding order, or negative if there was no such band.
 */
static inline int ff_sfdelta_can_remove_band(const SingleChannelElement *sce,
    const uint8_t *nextband, int prev_sf, int band)
{
    return prev_sf >= 0
        && sce->sf_idx[nextband[band]] >= (prev_sf - SCALE_MAX_DIFF)
        && sce->sf_idx[nextband[band]] <= (prev_sf + SCALE_MAX_DIFF);
}

/*
 * Checks whether the specified band's scalefactor could be replaced
 * with another one without violating SF delta encoding constraints.
 * prev_sf has to be the scalefactor of the previous nonzero, nonsepcial
 * band, in encoding order, or negative if there was no such band.
 */
static inline int ff_sfdelta_can_replace(const SingleChannelElement *sce,
    const uint8_t *nextband, int prev_sf, int new_sf, int band)
{
    return new_sf >= (prev_sf - SCALE_MAX_DIFF)
        && new_sf <= (prev_sf + SCALE_MAX_DIFF)
        && sce->sf_idx[nextband[band]] >= (new_sf - SCALE_MAX_DIFF)
        && sce->sf_idx[nextband[band]] <= (new_sf + SCALE_MAX_DIFF);
}

#define ERROR_IF(cond, ...) \
    if (cond) { \
        av_log(avctx, AV_LOG_ERROR, __VA_ARGS__); \
        return AVERROR(EINVAL); \
    }

#define WARN_IF(cond, ...) \
    if (cond) { \
        av_log(avctx, AV_LOG_WARNING, __VA_ARGS__); \
    }

#endif /* AVCODEC_AACENC_UTILS_H */
