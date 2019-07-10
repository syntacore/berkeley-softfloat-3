
/** @file

This C header file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
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


/** @file
@note  If SoftFloat is made available as a general library for programs to
use, it is strongly recommended that a platform-specific version of this
header, "softfloat.h", be created that folds in "softfloat_types.h" and that
eliminates all dependencies on compile-time macros.
*/


#ifndef SOFTFLOAT_FUNCTIONS_H_
#define SOFTFLOAT_FUNCTIONS_H_

#include "softfloat/types.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

    /**
    Routine to raise any or all of the software floating-point exception flags.
    */
    /**@{*/
    void softfloat_raiseFlags(uint8_t);
    void softfloat_clearFlags(void);
    uint8_t softfloat_getFlags(void);
    /**@}*/

    softfloat_round_mode softfloat_get_roundingMode(void);
    void softfloat_set_roundingMode(softfloat_round_mode);

    /**
    Integer-to-floating-point conversion routines.
    */
    /**@{*/
    float16_t ui32_to_f16(uint32_t);
    float32_t ui32_to_f32(uint32_t);
    float64_t ui32_to_f64(uint32_t);
    float16_t ui64_to_f16(uint64_t);
    float32_t ui64_to_f32(uint64_t);
    float64_t ui64_to_f64(uint64_t);
    float16_t i32_to_f16(int32_t);
    float32_t i32_to_f32(int32_t);
    float64_t i32_to_f64(int32_t);
    float16_t i64_to_f16(int64_t);
    float32_t i64_to_f32(int64_t);
    float64_t i64_to_f64(int64_t);
    /**@}*/

    /**
    16-bit (half-precision) floating-point operations.
    */
    /**@{*/
    uint32_t f16_to_ui32(float16_t, softfloat_round_mode, bool);
    uint64_t f16_to_ui64(float16_t, softfloat_round_mode, bool);
    int32_t f16_to_i32(float16_t, softfloat_round_mode, bool);
    int64_t f16_to_i64(float16_t, softfloat_round_mode, bool);
    uint32_t f16_to_ui32_r_minMag(float16_t, bool);
    uint64_t f16_to_ui64_r_minMag(float16_t, bool);
    int32_t f16_to_i32_r_minMag(float16_t, bool);
    int64_t f16_to_i64_r_minMag(float16_t, bool);
    float32_t f16_to_f32(float16_t);
    float64_t f16_to_f64(float16_t);
    float16_t f16_roundToInt(float16_t, uint8_t, bool);
    float16_t f16_add(float16_t, float16_t);
    float16_t f16_sub(float16_t, float16_t);
    float16_t f16_mul(float16_t, float16_t);
    float16_t f16_mulAdd(float16_t, float16_t, float16_t);
    float16_t f16_div(float16_t, float16_t);
    float16_t f16_rem(float16_t, float16_t);
    float16_t f16_sqrt(float16_t);
    bool f16_eq(float16_t, float16_t);
    bool f16_le(float16_t, float16_t);
    bool f16_lt(float16_t, float16_t);
    bool f16_eq_signaling(float16_t, float16_t);
    bool f16_le_quiet(float16_t, float16_t);
    bool f16_lt_quiet(float16_t, float16_t);
    bool f16_isSignalingNaN(float16_t);
    /**@}*/

    /**
    32-bit (single-precision) floating-point operations.
    */
    /**@{*/
    uint32_t f32_to_ui32(float32_t, softfloat_round_mode, bool);
    uint64_t f32_to_ui64(float32_t, softfloat_round_mode, bool);
    int32_t f32_to_i32(float32_t, softfloat_round_mode, bool);
    int64_t f32_to_i64(float32_t, softfloat_round_mode, bool);
    uint32_t f32_to_ui32_r_minMag(float32_t, bool);
    uint64_t f32_to_ui64_r_minMag(float32_t, bool);
    int32_t f32_to_i32_r_minMag(float32_t, bool);
    int64_t f32_to_i64_r_minMag(float32_t, bool);
    float16_t f32_to_f16(float32_t);
    float64_t f32_to_f64(float32_t);
    float32_t f32_roundToInt(float32_t, uint8_t, bool);
    float32_t f32_add(float32_t, float32_t);
    float32_t f32_sub(float32_t, float32_t);
    float32_t f32_mul(float32_t, float32_t);
    float32_t f32_mulAdd(float32_t, float32_t, float32_t);
    float32_t f32_div(float32_t, float32_t);
    float32_t f32_rem(float32_t, float32_t);
    float32_t f32_sqrt(float32_t);
    bool f32_eq(float32_t, float32_t);
    bool f32_le(float32_t, float32_t);
    bool f32_lt(float32_t, float32_t);
    bool f32_eq_signaling(float32_t, float32_t);
    bool f32_le_quiet(float32_t, float32_t);
    bool f32_lt_quiet(float32_t, float32_t);
    bool f32_isSignalingNaN(float32_t);
    /**@}*/

    /** 64-bit (double-precision) floating-point operations. */
    /**@{*/
    uint32_t f64_to_ui32(float64_t, softfloat_round_mode, bool);
    uint64_t f64_to_ui64(float64_t, softfloat_round_mode, bool);
    int32_t f64_to_i32(float64_t, softfloat_round_mode, bool);
    int64_t f64_to_i64(float64_t, softfloat_round_mode, bool);
    uint32_t f64_to_ui32_r_minMag(float64_t, bool);
    uint64_t f64_to_ui64_r_minMag(float64_t, bool);
    int32_t f64_to_i32_r_minMag(float64_t, bool);
    int64_t f64_to_i64_r_minMag(float64_t, bool);
    float16_t f64_to_f16(float64_t);
    float32_t f64_to_f32(float64_t);
    void f64_to_extF80M(float64_t, extFloat80_t*);
    float64_t f64_roundToInt(float64_t, uint8_t, bool);
    float64_t f64_add(float64_t, float64_t);
    float64_t f64_sub(float64_t, float64_t);
    float64_t f64_mul(float64_t, float64_t);
    float64_t f64_mulAdd(float64_t, float64_t, float64_t);
    float64_t f64_div(float64_t, float64_t);
    float64_t f64_rem(float64_t, float64_t);
    float64_t f64_sqrt(float64_t);
    bool f64_eq(float64_t, float64_t);
    bool f64_le(float64_t, float64_t);
    bool f64_lt(float64_t, float64_t);
    bool f64_eq_signaling(float64_t, float64_t);
    bool f64_le_quiet(float64_t, float64_t);
    bool f64_lt_quiet(float64_t, float64_t);
    bool f64_isSignalingNaN(float64_t);
    /**@}*/

    bool f128M_isSignalingNaN(const float128_t*);

    void extF80M_add(extFloat80_t const*, extFloat80_t const*, extFloat80_t*);
    void f128M_add(float128_t const*, float128_t const*, float128_t*);

    void extF80M_sub(extFloat80_t const*, extFloat80_t const*, extFloat80_t*);
    void f128M_sub(float128_t const*, float128_t const*, float128_t*);

    void f128M_mulAdd(float128_t const*, float128_t const*, float128_t const*, float128_t*);

#ifndef SOFTFLOAT_FAST_INT64

    void i32_to_extF80M(int32_t, extFloat80_t*);
    void i32_to_f128M(int32_t, float128_t*);
    void f16_to_f128M(float16_t, float128_t*);
    void f64_to_f128M(float64_t, float128_t*);
    void f32_to_extF80M(float32_t, extFloat80_t*);
    void f16_to_extF80M(float16_t, extFloat80_t*);
    void f32_to_f128M(float32_t, float128_t*);
    void ui64_to_extF80M(uint64_t, extFloat80_t*);
    void i64_to_extF80M(int64_t, extFloat80_t*);
    void ui64_to_f128M(uint64_t, float128_t*);
    void i64_to_f128M(int64_t, float128_t*);
    void ui32_to_f128M(uint32_t, float128_t*);
    void ui32_to_extF80M(uint32_t, extFloat80_t*);
    uint32_t extF80M_to_ui32(extFloat80_t const*, softfloat_round_mode, bool);
    uint64_t extF80M_to_ui64(extFloat80_t const*, softfloat_round_mode, bool);
    int32_t extF80M_to_i32(extFloat80_t const*, softfloat_round_mode, bool);
    int64_t extF80M_to_i64(extFloat80_t const*, softfloat_round_mode, bool);
    uint32_t extF80M_to_ui32_r_minMag(extFloat80_t const*, bool);
    uint64_t extF80M_to_ui64_r_minMag(extFloat80_t const*, bool);
    int32_t extF80M_to_i32_r_minMag(extFloat80_t const*, bool);
    int64_t extF80M_to_i64_r_minMag(extFloat80_t const*, bool);
    float16_t extF80M_to_f16(extFloat80_t const*);
    float32_t extF80M_to_f32(extFloat80_t const*);
    float64_t extF80M_to_f64(extFloat80_t const*);
    void extF80M_to_f128M(extFloat80_t const*, float128_t*);
    void extF80M_roundToInt(extFloat80_t const*, uint8_t, bool, extFloat80_t*);
    void extF80M_mul(extFloat80_t const*, extFloat80_t const*, extFloat80_t*);
    void extF80M_div(extFloat80_t const*, extFloat80_t const*, extFloat80_t*);
    void extF80M_rem(extFloat80_t const*, extFloat80_t const*, extFloat80_t*);
    void extF80M_sqrt(extFloat80_t const*, extFloat80_t*);
    bool extF80M_eq(extFloat80_t const*, extFloat80_t const*);
    bool extF80M_le(extFloat80_t const*, extFloat80_t const*);
    bool extF80M_lt(extFloat80_t const*, extFloat80_t const*);
    bool extF80M_eq_signaling(const extFloat80_t*, const extFloat80_t*);
    bool extF80M_le_quiet(const extFloat80_t*, const extFloat80_t*);
    bool extF80M_lt_quiet(const extFloat80_t*, const extFloat80_t*);

    uint32_t f128M_to_ui32(float128_t const*, softfloat_round_mode, bool);
    uint64_t f128M_to_ui64(float128_t const*, softfloat_round_mode, bool);
    int32_t f128M_to_i32(float128_t const*, softfloat_round_mode, bool);
    int64_t f128M_to_i64(float128_t const*, softfloat_round_mode, bool);
    uint32_t f128M_to_ui32_r_minMag(const float128_t*, bool);
    uint64_t f128M_to_ui64_r_minMag(const float128_t*, bool);
    int32_t f128M_to_i32_r_minMag(const float128_t*, bool);
    int64_t f128M_to_i64_r_minMag(const float128_t*, bool);
    float16_t f128M_to_f16(const float128_t*);
    float32_t f128M_to_f32(const float128_t*);
    float64_t f128M_to_f64(const float128_t*);
    void f128M_to_extF80M(const float128_t*, extFloat80_t*);
    void f128M_roundToInt(const float128_t*, uint8_t, bool, float128_t*);
    void f128M_mul(const float128_t*, const float128_t*, float128_t*);
    void f128M_div(const float128_t*, const float128_t*, float128_t*);
    void f128M_rem(const float128_t*, const float128_t*, float128_t*);
    void f128M_sqrt(const float128_t*, float128_t*);
    bool f128M_eq(const float128_t*, const float128_t*);
    bool f128M_le(const float128_t*, const float128_t*);
    bool f128M_lt(const float128_t*, const float128_t*);
    bool f128M_eq_signaling(const float128_t*, const float128_t*);
    bool f128M_le_quiet(const float128_t*, const float128_t*);
    bool f128M_lt_quiet(const float128_t*, const float128_t*);

#else

    extFloat80_t ui32_to_extF80(uint32_t);
    extFloat80_t ui64_to_extF80(uint64_t);
    extFloat80_t i32_to_extF80(int32_t);
    extFloat80_t i64_to_extF80(int64_t);
    extFloat80_t f16_to_extF80(float16_t);
    extFloat80_t f32_to_extF80(float32_t);
    extFloat80_t f64_to_extF80(float64_t);

    uint32_t extF80_to_ui32(extFloat80_t, softfloat_round_mode, bool);
    uint64_t extF80_to_ui64(extFloat80_t, softfloat_round_mode, bool);
    int32_t extF80_to_i32(extFloat80_t, softfloat_round_mode, bool);
    int64_t extF80_to_i64(extFloat80_t, softfloat_round_mode, bool);
    uint32_t extF80_to_ui32_r_minMag(extFloat80_t, bool);
    uint64_t extF80_to_ui64_r_minMag(extFloat80_t, bool);
    int32_t extF80_to_i32_r_minMag(extFloat80_t, bool);
    int64_t extF80_to_i64_r_minMag(extFloat80_t, bool);

    float16_t extF80_to_f16(extFloat80_t);
    float32_t extF80_to_f32(extFloat80_t);
    float64_t extF80_to_f64(extFloat80_t);
    float128_t extF80_to_f128(extFloat80_t);

    extFloat80_t extF80_roundToInt(extFloat80_t, softfloat_round_mode, bool);
    extFloat80_t extF80_add(extFloat80_t, extFloat80_t);
    extFloat80_t extF80_sub(extFloat80_t, extFloat80_t);
    extFloat80_t extF80_mul(extFloat80_t, extFloat80_t);
    extFloat80_t extF80_div(extFloat80_t, extFloat80_t);
    extFloat80_t extF80_rem(extFloat80_t, extFloat80_t);
    extFloat80_t extF80_sqrt(extFloat80_t);

    bool extF80_eq(extFloat80_t, extFloat80_t);
    bool extF80_le(extFloat80_t, extFloat80_t);
    bool extF80_lt(extFloat80_t, extFloat80_t);
    bool extF80_eq_signaling(extFloat80_t, extFloat80_t);
    bool extF80_le_quiet(extFloat80_t, extFloat80_t);
    bool extF80_lt_quiet(extFloat80_t, extFloat80_t);
    bool extF80_isSignalingNaN(extFloat80_t);

    float128_t ui32_to_f128(uint32_t);
    float128_t ui64_to_f128(uint64_t);
    float128_t i32_to_f128(int32_t);
    float128_t i64_to_f128(int64_t);
    float128_t f16_to_f128(float16_t);
    float128_t f32_to_f128(float32_t);
    float128_t f64_to_f128(float64_t);

    uint32_t f128_to_ui32(float128_t, softfloat_round_mode, bool);
    uint64_t f128_to_ui64(float128_t, softfloat_round_mode, bool);
    int32_t f128_to_i32(float128_t, softfloat_round_mode, bool);
    int64_t f128_to_i64(float128_t, softfloat_round_mode, bool);

    uint32_t f128_to_ui32_r_minMag(float128_t, bool);
    uint64_t f128_to_ui64_r_minMag(float128_t, bool);
    int32_t f128_to_i32_r_minMag(float128_t, bool);
    int64_t f128_to_i64_r_minMag(float128_t, bool);

    float16_t f128_to_f16(float128_t);
    float32_t f128_to_f32(float128_t);
    float64_t f128_to_f64(float128_t);
    extFloat80_t f128_to_extF80(float128_t);

    float128_t f128_roundToInt(float128_t, softfloat_round_mode, bool);
    float128_t f128_add(float128_t, float128_t);
    float128_t f128_sub(float128_t, float128_t);
    float128_t f128_mul(float128_t, float128_t);
    float128_t f128_mulAdd(float128_t, float128_t, float128_t);
    float128_t f128_div(float128_t, float128_t);
    float128_t f128_rem(float128_t, float128_t);
    float128_t f128_sqrt(float128_t);

    bool f128_eq(float128_t, float128_t);
    bool f128_le(float128_t, float128_t);
    bool f128_lt(float128_t, float128_t);
    bool f128_eq_signaling(float128_t, float128_t);
    bool f128_le_quiet(float128_t, float128_t);
    bool f128_lt_quiet(float128_t, float128_t);
    bool f128_isSignalingNaN(float128_t);

    static inline void
        i32_to_extF80M(int32_t a,
                       extFloat80_t* zPtr)
    {
        *zPtr = i32_to_extF80(a);
    }

    static inline void
        i32_to_f128M(int32_t a, float128_t* zPtr)
    {
        *zPtr = i32_to_f128(a);
    }

    static inline void
        f16_to_f128M(float16_t a,
                     float128_t* zPtr)
    {
        *zPtr = f16_to_f128(a);
    }

    static inline void
        f64_to_f128M(float64_t a,
                     float128_t* zPtr)
    {
        *zPtr = f64_to_f128(a);
    }

    static inline void
        f32_to_extF80M(float32_t a, extFloat80_t* zPtr)
    {
        *zPtr = f32_to_extF80(a);
    }

    static inline void
        f16_to_extF80M(float16_t a, extFloat80_t* zPtr)
    {
        *zPtr = f16_to_extF80(a);
    }

    static inline void
        f32_to_f128M(float32_t a, float128_t* zPtr)
    {
        *zPtr = f32_to_f128(a);
    }

    static inline void
        ui64_to_extF80M(uint64_t a,
                        extFloat80_t* const zPtr)
    {
        *zPtr = ui64_to_extF80(a);
    }

    static inline void
        i64_to_extF80M(int64_t a, extFloat80_t* zPtr)
    {
        *zPtr = i64_to_extF80(a);
    }

    static inline void
        ui64_to_f128M(uint64_t a,
                      float128_t* const zPtr)
    {
        *zPtr = ui64_to_f128(a);
    }

    static inline void
        i64_to_f128M(int64_t a,
                     float128_t* const zPtr)
    {
        *zPtr = i64_to_f128(a);
    }

    static inline void
        ui32_to_extF80M(uint32_t a, extFloat80_t* zPtr)
    {
        *zPtr = ui32_to_extF80(a);
    }

    static inline void
        ui32_to_f128M(uint32_t a, float128_t* zPtr)
    {
        *zPtr = ui32_to_f128(a);
    }

    static inline uint32_t
        extF80M_to_ui32(extFloat80_t const* const aPtr,
                        softfloat_round_mode const roundingMode,
                        bool exact)
    {
        return extF80_to_ui32(*aPtr, roundingMode, exact);
    }

    static inline uint64_t
        extF80M_to_ui64(const extFloat80_t* aPtr,
                        softfloat_round_mode const roundingMode,
                        bool exact)
    {
        return extF80_to_ui64(*aPtr, roundingMode, exact);
    }

    static inline int32_t
        extF80M_to_i32(extFloat80_t const* const aPtr,
                       softfloat_round_mode const roundingMode,
                       bool const exact)
    {
        return extF80_to_i32(*aPtr, roundingMode, exact);
    }

    static inline int64_t
        extF80M_to_i64(extFloat80_t const* const aPtr,
                       softfloat_round_mode const roundingMode,
                       bool const exact)
    {
        return extF80_to_i64(*aPtr, roundingMode, exact);
    }

    static inline uint32_t
        extF80M_to_ui32_r_minMag(extFloat80_t const* const aPtr,
                                 bool const exact)
    {
        return extF80_to_ui32_r_minMag(*aPtr, exact);
    }

    static inline uint64_t
        extF80M_to_ui64_r_minMag(extFloat80_t const* const aPtr,
                                 bool const exact)
    {
        return extF80_to_ui64_r_minMag(*aPtr, exact);
    }

    static inline int32_t
        extF80M_to_i32_r_minMag(extFloat80_t const* const aPtr,
                                bool const exact)
    {
        return extF80_to_i32_r_minMag(*aPtr, exact);
    }

    static inline int64_t
        extF80M_to_i64_r_minMag(extFloat80_t const* const aPtr,
                                bool const exact)
    {
        return extF80_to_i64_r_minMag(*aPtr, exact);
    }

    static inline float16_t
        extF80M_to_f16(extFloat80_t const* const aPtr)
    {
        return extF80_to_f16(*aPtr);
    }

    static inline float32_t
        extF80M_to_f32(extFloat80_t const* const aPtr)
    {
        return extF80_to_f32(*aPtr);
    }

    static inline float64_t
        extF80M_to_f64(const extFloat80_t* aPtr)
    {
        return extF80_to_f64(*aPtr);
    }

    static inline void
        extF80M_to_f128M(const extFloat80_t* aPtr, float128_t* zPtr)
    {
        *zPtr = extF80_to_f128(*aPtr);
    }

    static inline void
        extF80M_roundToInt(extFloat80_t const* aPtr,
                           softfloat_round_mode const roundingMode,
                           bool const exact,
                           extFloat80_t* zPtr)
    {
        *zPtr = extF80_roundToInt(*aPtr, roundingMode, exact);
    }

    static inline void
        extF80M_mul(extFloat80_t const* const aPtr,
                    extFloat80_t const* const  bPtr,
                    extFloat80_t* const zPtr)
    {
        *zPtr = extF80_mul(*aPtr, *bPtr);
    }

    static inline void
        extF80M_div(extFloat80_t const* const aPtr,
                    extFloat80_t const* const bPtr,
                    extFloat80_t* const zPtr)
    {
        *zPtr = extF80_div(*aPtr, *bPtr);
    }

    static inline void
        extF80M_rem(extFloat80_t const* const aPtr,
                    extFloat80_t const* const bPtr,
                    extFloat80_t* const zPtr)
    {
        *zPtr = extF80_rem(*aPtr, *bPtr);
    }

    static inline void
        extF80M_sqrt(extFloat80_t const* const aPtr,
                     extFloat80_t* const zPtr)
    {
        *zPtr = extF80_sqrt(*aPtr);
    }

    static inline bool
        extF80M_eq(extFloat80_t const* const aPtr,
                   extFloat80_t const* const bPtr)
    {
        return extF80_eq(*aPtr, *bPtr);
    }

    static inline bool
        extF80M_le(extFloat80_t const* const aPtr,
                   extFloat80_t const* const bPtr)
    {
        return extF80_le(*aPtr, *bPtr);
    }

    static inline bool
        extF80M_lt(extFloat80_t const* const aPtr,
                   extFloat80_t const* const bPtr)
    {
        return extF80_lt(*aPtr, *bPtr);
    }

    static inline bool
        extF80M_eq_signaling(extFloat80_t const* const aPtr,
                             extFloat80_t const* const bPtr)
    {
        return extF80_eq_signaling(*aPtr, *bPtr);
    }

    static inline bool
        extF80M_le_quiet(extFloat80_t const* const aPtr,
                         extFloat80_t const* const bPtr)
    {
        return extF80_le_quiet(*aPtr, *bPtr);
    }

    static inline bool
        extF80M_lt_quiet(extFloat80_t const* const aPtr,
                         extFloat80_t const* const bPtr)
    {
        return extF80_lt_quiet(*aPtr, *bPtr);
    }


    static inline uint32_t
        f128M_to_ui32(float128_t const* const aPtr,
                      softfloat_round_mode const roundingMode,
                      bool const exact)
    {
        return f128_to_ui32(*aPtr, roundingMode, exact);
    }

    static inline uint64_t
        f128M_to_ui64(float128_t const* const aPtr,
                      softfloat_round_mode const roundingMode,
                      bool const exact)
    {
        return f128_to_ui64(*aPtr, roundingMode, exact);
    }

    static inline int32_t
        f128M_to_i32(float128_t const* const aPtr,
                     softfloat_round_mode const roundingMode,
                     bool const exact)
    {
        return f128_to_i32(*aPtr, roundingMode, exact);
    }

    static inline int64_t
        f128M_to_i64(float128_t const* const aPtr,
                     softfloat_round_mode const roundingMode,
                     bool const exact)
    {
        return f128_to_i64(*aPtr, roundingMode, exact);
    }

    static inline uint32_t
        f128M_to_ui32_r_minMag(float128_t const* aPtr,
                               bool exact)
    {
        return f128_to_ui32_r_minMag(*aPtr, exact);
    }

    static inline uint64_t
        f128M_to_ui64_r_minMag(const float128_t* aPtr, bool exact)
    {
        return f128_to_ui64_r_minMag(*aPtr, exact);
    }

    static inline int32_t
        f128M_to_i32_r_minMag(const float128_t* aPtr,
                              bool exact)
    {
        return f128_to_i32_r_minMag(*aPtr, exact);
    }

    static inline int64_t
        f128M_to_i64_r_minMag(const float128_t* aPtr,
                              bool exact)
    {
        return f128_to_i64_r_minMag(*aPtr, exact);
    }

    static inline float16_t
        f128M_to_f16(const float128_t* aPtr)
    {
        return f128_to_f16(*aPtr);
    }

    static inline float32_t
        f128M_to_f32(const float128_t* aPtr)
    {
        return f128_to_f32(*aPtr);
    }

    static inline float64_t
        f128M_to_f64(const float128_t* aPtr)
    {
        return f128_to_f64(*aPtr);
    }

    static inline void
        f128M_to_extF80M(const float128_t* aPtr,
                         extFloat80_t* zPtr)
    {
        *zPtr = f128_to_extF80(*aPtr);
    }

    static inline void
        f128M_roundToInt(float128_t const* aPtr,
                         softfloat_round_mode const roundingMode,
                         bool const exact,
                         float128_t* zPtr)
    {
        *zPtr = f128_roundToInt(*aPtr, roundingMode, exact);
    }

    static inline void
        f128M_mul(const float128_t* aPtr,
                  const float128_t* bPtr,
                  float128_t* zPtr)
    {
        *zPtr = f128_mul(*aPtr, *bPtr);
    }

    static inline void
        f128M_div(const float128_t* aPtr, const float128_t* bPtr, float128_t* zPtr)
    {
        *zPtr = f128_div(*aPtr, *bPtr);
    }

    static inline void
        f128M_rem(const float128_t* aPtr,
                  const float128_t* bPtr,
                  float128_t* zPtr)
    {
        *zPtr = f128_rem(*aPtr, *bPtr);
    }

    static inline void
        f128M_sqrt(const float128_t* aPtr,
                   float128_t* zPtr)
    {
        *zPtr = f128_sqrt(*aPtr);
    }

    static inline bool
        f128M_eq(const float128_t* aPtr, const float128_t* bPtr)
    {
        return f128_eq(*aPtr, *bPtr);
    }

    static inline bool
        f128M_le(const float128_t* aPtr,
                 const float128_t* bPtr)
    {
        return f128_le(*aPtr, *bPtr);
    }

    static inline bool
        f128M_lt(const float128_t* aPtr,
                 const float128_t* bPtr)
    {
        return f128_lt(*aPtr, *bPtr);
    }

    static inline bool
        f128M_eq_signaling(const float128_t* aPtr, const float128_t* bPtr)
    {
        return f128_eq_signaling(*aPtr, *bPtr);
    }

    static inline bool
        f128M_le_quiet(const float128_t* aPtr,
                       const float128_t* bPtr)
    {
        return f128_le_quiet(*aPtr, *bPtr);
    }

    static inline bool
        f128M_lt_quiet(const float128_t* aPtr,
                       const float128_t* bPtr)
    {
        return f128_lt_quiet(*aPtr, *bPtr);
    }

#endif

    inline bool
        extF80M_isSignalingNaN(extFloat80_t const* const aPtr)
    {
        return
            UINT16_C(0x7FFF) == (UINT16_C(0x7FFF) & aPtr->signExp) &&
            0 == (UINT64_C(0x4000000000000000) & aPtr->signif) &&
            0 != (UINT64_C(0x3FFFFFFFFFFFFFFF) & aPtr->signif);
    }

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  /* SOFTFLOAT_FUNCTIONS_H_ */

