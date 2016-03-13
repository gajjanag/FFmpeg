/*
 * 64 bit random number generator written by
 * Written in 2015 by Sebastiano Vigna (vigna@acm.org)
 * under public domain:
 * https://creativecommons.org/publicdomain/zero/1.0/
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

#ifndef AVUTIL_RAND_H
#define AVUTIL_RAND_H
#include <stdint.h>

typedef struct AVRAND64 {
    uint64_t state[2];
} AVRAND64;

/**
 * Initialize the 64 bit random number generator.
 *
 * @param rng AVRAND64 structure that is initialized
 * @param seed 64 bit seed
 */
void av_rand64_init(AVRAND64 *rng, uint64_t seed);

/**
 * Get the next 64 bit random number from the rng.
 *
 * @param rng AVRAND64 structure holding the state of the rng
 * @return 64 bit random number
 */
static inline uint64_t av_rand64_get(AVRAND64 *rng){
    uint64_t x = rng->state[0];
    uint64_t const y = rng->state[1];
    rng->state[0] = y;
    x ^= x << 23; // a
    rng->state[1] = x ^ y ^ (x >> 17) ^ (y >> 26); // b, c
    return rng->state[1] + y;
}

#endif /* AVUTIL_RAND_H */
