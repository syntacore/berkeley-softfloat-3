/** @file

This C header file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015 The Regents of the University of
California.  All rights reserved.
*/
/*
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions, and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions, and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 3. Neither the name of the University nor the names of its contributors may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS", AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, ARE
DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef SOFTFLOAT_TYPES_H
#define SOFTFLOAT_TYPES_H

#include <stdint.h>

#ifndef THREAD_LOCAL
#define THREAD_LOCAL
#endif

#ifdef __cplusplus
extern "C" {
#endif

    /**
    Types used to pass 16-bit, 32-bit, 64-bit, and 128-bit floating-point
    arguments and results to/from functions.  These types must be exactly
    16 bits, 32 bits, 64 bits, and 128 bits in size, respectively.  Where a
    platform has "native" support for IEEE-Standard floating-point formats,
    the types below may, if desired, be defined as aliases for the native types
    (typically `float' and `double', and possibly `long double').
    */
    struct float16_type
    {
        uint16_t v;
    };
    typedef struct float16_type float16_t;

    struct float32_type
    {
        uint32_t v;
    };

    typedef struct float32_type float32_t;

    struct float64_type
    {
        uint64_t v;
    };
    
    typedef struct float64_type float64_t;

    struct float128_type
    {
#ifdef BIG_ENDIAN
        uint64_t v64;
        uint64_t v0;
#else
        uint64_t v0;
        uint64_t v64;
#endif
    };

    typedef struct float128_type float128_t;

    /**
    The format of an 80-bit extended floating-point number in memory.  This
    structure must contain a 16-bit field named `signExp' and a 64-bit field
    named `signif'.
    */
    /** @bug macro value dependent interface */
    struct extFloat80M
    {
#ifdef BIG_ENDIAN
        uint16_t signExp;
        uint64_t signif;
#else
        uint64_t signif;
        uint16_t signExp;
#endif
    };

    /**
    The type used to pass 80-bit extended floating-point arguments and
    results to/from functions.  This type must have size identical to
    `struct extFloat80M'.  Type `extFloat80_t' can be defined as an alias for
    `struct extFloat80M'.  Alternatively, if a platform has "native" support
    for IEEE-Standard 80-bit extended floating-point, it may be possible,
    if desired, to define `extFloat80_t' as an alias for the native type
    (presumably either `long double' or a nonstandard compiler-intrinsic type).
    In that case, the `signif' and `signExp' fields of `struct extFloat80M'
    must align exactly with the locations in memory of the sign, exponent, and
    significand of the native type.
    */
    typedef struct extFloat80M extFloat80_t;

    /**
    Software floating-point exception flags.
    */
    enum softfloat_flags
    {
        softfloat_flag_inexact = 1 << 0,
        softfloat_flag_underflow = 1 << 1,
        softfloat_flag_overflow = 1 << 2,
        softfloat_flag_infinite = 1 << 3,
        softfloat_flag_invalid = 1 << 4
    };

    /**
    Software floating-point underflow tininess-detection mode.
    */
    typedef enum softfloat_tininess
    {
        softfloat_tininess_beforeRounding = 0,
        softfloat_tininess_afterRounding = 1
    } softfloat_tininess;

    /**
    Software floating-point rounding mode.
    */
    typedef enum softfloat_round_mode
    {
        softfloat_round_near_even = 0,
        softfloat_round_minMag = 1,
        softfloat_round_min = 2,
        softfloat_round_max = 3,
        softfloat_round_near_maxMag = 4
    } softfloat_round_mode;

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // SOFTFLOAT_TYPES_H
