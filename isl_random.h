#ifndef ISL_INCLUDE_RANDOM_H_
#define ISL_INCLUDE_RANDOM_H_

/* isl_random - v0.1 public domain pseudo random number generator  
                        no warranty implied; use at your own risk
   
   Do this:
       #define ISL_RANDOM_IMPLEMENTATION
   before you include this file in *one* C or C++ file to create the implementation.

   To static link also add:
       #define ISL_RANDOM_STATIC

   QUICK NOTES:
	    This is just a simple wrapper around Xorshiro256**(XOR, shift, rotate) library
			taken from https://prng.di.unimi.it/xoshiro256starstar.c which is licensed
			under CC0 license (see end of file). The state is not static but passed as an
			argument to the functions. Also there are some renamings to reduce namespace
			pollutions.

	 USAGE:
	    uint64_t state[ISL_RANDOM_STATE_SIZE];
	    isl_random_init(&state, 0xDEADBEEF);    // Use builtin Splitmix64 to init state
	    uint64_t raw = isl_random_next(&state);
	    int random_int = isl_random_int(&state, 0, 1000); // Generate random int [0-1000)
	    double random_double = isl_random_double(&state); // Generate random double [0.0-1.0)
	    printf("%d %5.5f\n", random_int, random_double);  // Should print 792 0.33190

   author: Ilya Kolbin (iskolbin@gmail.com)
   url: https://github.com/iskolbin/isl_random 
   git: git@github.com:iskolbin/isl_random

   LICENSE:
     See end of file for license information.
*/

/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>

/* This is xoshiro256** 1.0, one of our all-purpose, rock-solid
   generators. It has excellent (sub-ns) speed, a state (256 bits) that is
   large enough for any parallel application, and it passes all tests we
   are aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

#ifdef ISL_RANDOM_STATIC
#define ISL_RANDOM_DEF static
#else
#define ISL_RANDOM_DEF extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define ISL_RANDOM_STATE_SIZE 4

ISL_RANDOM_DEF void isl_random_init(uint64_t *state, uint64_t seed);
ISL_RANDOM_DEF double isl_random_double(uint64_t *state);
ISL_RANDOM_DEF int isl_random_int(uint64_t *state, int from, int to);

ISL_RANDOM_DEF uint64_t isl_random_next(uint64_t *state);
ISL_RANDOM_DEF void isl_random_jump(uint64_t *state);
ISL_RANDOM_DEF void isl_random_long_jump(uint64_t *state);

#ifdef __cplusplus
}
#endif

#ifdef ISL_RANDOM_IMPLEMENTATION
#ifndef ISL_RANDOM_IMPLEMENTATION_ONCE
#define ISL_RANDOM_IMPLEMENTATION_ONCE

static inline uint64_t isl_random__rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

ISL_RANDOM_DEF void isl_random_init(uint64_t *state, uint64_t seed) {
	for (int i = 0; i < ISL_RANDOM_STATE_SIZE; i++) {  /* Splitmix64 taken from Rosetta code */
		seed += 0x9e3779b97f4a7c15;                /* increment the state variable */
		uint64_t z = seed;                         /* copy the state to a working variable */
		z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;  /* xor the variable with the variable right bit shifted 30 then multiply by a constant */
		z = (z ^ (z >> 27)) * 0x94d049bb133111eb;  /* xor the variable with the variable right bit shifted 27 then multiply by a constant */
		state[i] = z ^ (z >> 31);                  /* return the variable xored with itself right bit shifted 31 */
	}
}

ISL_RANDOM_DEF uint64_t isl_random_next(uint64_t *state) {
	const uint64_t result = isl_random__rotl(state[1] * 5, 7) * 9;

	const uint64_t t = state[1] << 17;

	state[2] ^= state[0];
	state[3] ^= state[1];
	state[1] ^= state[2];
	state[0] ^= state[3];

	state[2] ^= t;

	state[3] = isl_random__rotl(state[3], 45);

	return result;
}

ISL_RANDOM_DEF double isl_random_double(uint64_t *state) {
	double y = (double) isl_random_next(state);
	return y / (0x8000000000000000U * 2.0);
}

ISL_RANDOM_DEF int isl_random_int(uint64_t *state, int from, int to) {
	if (from == to) return from;
	int d = (from > to) ? from - to : to - from;
	uint64_t v = isl_random_next(state);
	return (int) (v % d) + from;
}

/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

ISL_RANDOM_DEF void isl_random_jump(uint64_t *state) {
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= state[0];
				s1 ^= state[1];
				s2 ^= state[2];
				s3 ^= state[3];
			}
			isl_random_next(state);
		}

	state[0] = s0;
	state[1] = s1;
	state[2] = s2;
	state[3] = s3;
}

/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations. */

ISL_RANDOM_DEF void isl_random_long_jump(uint64_t *state) {
	static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= state[0];
				s1 ^= state[1];
				s2 ^= state[2];
				s3 ^= state[3];
			}
			isl_random_next(state);
		}

	state[0] = s0;
	state[1] = s1;
	state[2] = s2;
	state[3] = s3;
}
#endif
#endif
#endif
/*
------------------------------------------------------------------------------
Original xorshiro256** is licenesed under CC0 license
------------------------------------------------------------------------------
CC0 License

No Copyright

The person who associated a work with this deed has dedicated the work to the
public domain by waiving all of his or her rights to the work worldwide under
copyright law, including all related and neighboring rights, to the extent
allowed by law.

You can copy, modify, distribute and perform the work, even for commercial
purposes, all without asking permission. See Other Information below.

This license is acceptable for Free Cultural Works.

Other Information

* In no way are the patent or trademark rights of any person affected by CC0,
nor are the rights that other persons may have in the work or in how the work
is used, such as publicity or privacy rights.
* Unless expressly stated otherwise, the person who associated a work with this
deed makes no warranties about the work, and disclaims liability for all uses
of the work, to the fullest extent permitted by applicable law.
* When using or citing the work, you should not imply endorsement by the author
or the affirmer.
*/

/*
------------------------------------------------------------------------------
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - MIT License
Copyright (c) 2022 Ilya Kolbin
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------
*/
