/** @file

This C header file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

@copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All rights reserved.

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
*/
/*
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

#ifndef SOFTFLOAT_INTERNALS_H_
#define SOFTFLOAT_INTERNALS_H_

#include "primitives/functions.hpp"
#include "softfloat/types.h"

#include <cassert>

#define signF16UI( a ) ((bool) ((uint16_t) (a)>>15))
#define expF16UI( a ) ((int8_t) ((a)>>10) & 0x1F)
#define fracF16UI( a ) ((a) & 0x03FF)
#define packToF16UI( sign, exp, sig ) (((uint16_t) (sign)<<15) + ((uint16_t) (exp)<<10) + (sig))

#define isNaNF16UI( a ) (((~(a) & 0x7C00) == 0) && ((a) & 0x03FF))

static inline uint16_t
f_as_u_16(float16_t v)
{
    return *(uint16_t const*)&v;
}
static inline float16_t
u_as_f_16(uint16_t v)
{
    return *(float16_t const*)&v;
}
static inline uint32_t
f_as_u_32(float32_t v)
{
    return *(uint32_t const*)&v;
}
static inline float32_t
u_as_f_32(uint32_t v)
{
    return *(float32_t const*)&v;
}

static inline uint64_t
f_as_u_64(float64_t v)
{
    return *(uint64_t const*)&v;
}
static inline float64_t
u_as_f_64(uint64_t v)
{
    return *(float64_t const*)&v;
}

static inline bool
signF32UI(uint32_t a)
{
    return (int32_t)a < 0;
}

/** @bug return signed 16-bits value instead unsigned 8-bits value */
static inline int16_t
expF32UI(uint32_t a)
{
    return (int16_t)((a >> 23) & ~(~UINT32_C(0) << 8));
}

static inline uint32_t
fracF32UI(uint32_t a)
{
    return a & ~(~UINT32_C(0) << 23);
}

static inline uint32_t
packToF32UI(bool sign, int16_t exp, uint32_t sgnf)
{
#if 0
    assert(0 <= exp && exp <= 255);
    assert(0 == (sgnf & (~UINT32_C(0) << 23)));
#endif
    uint32_t const u_res = ((uint32_t)exp << 23) + sgnf;
    assert(0 == (u_res & (UINT32_C(1) << 31)));
    return ((uint32_t)!!sign << 31) | u_res;
}

static inline uint32_t
signed_zero_F32UI(bool sign)
{
    return packToF32UI(sign, 0, 0);
}

static inline float32_t
signed_zero_F32(bool sign)
{
    return u_as_f_32(signed_zero_F32UI(sign));
}

static inline uint32_t
signed_inf_F32UI(bool sign)
{
    return packToF32UI(sign, 0xFF, 0);
}

static inline float32_t
signed_inf_F32(bool sign)
{
    return u_as_f_32(signed_inf_F32UI(sign));
}

static inline bool
isNaNF32UI(uint32_t a)
{
    return 255 == expF32UI(a) && 0 != fracF32UI(a);
}

static inline bool
isInf32UI(uint32_t a)
{
    return 255 == expF32UI(a) && 0 == fracF32UI(a);
}

static inline bool
isZero32UI(uint32_t a)
{
    return 0 == expF32UI(a) && 0 == fracF32UI(a);
}

static inline bool
isSubnormal32UI(uint32_t a)
{
    return 0 == expF32UI(a) && 0 != fracF32UI(a);
}

static inline bool
signF64UI(uint64_t a)
{
    return (int64_t)a < 0;
}
static int16_t
expF64UI(uint64_t a)
{
    return (int16_t)((a >> 52) & ~(~UINT32_C(0) << 11));
}
static inline uint64_t
fracF64UI(uint64_t a)
{
    return a & ~(~UINT64_C(0) << 52);
}
static inline uint64_t
packToF64UI(bool sign, int16_t exp, uint64_t sgnf)
{
#if 0
    assert(0 <= exp && exp <= 2047);
    assert(0 == (sgnf & (~UINT32_C(0) << 52)));
#endif
    uint64_t const u_res = ((uint64_t)exp << 52) + sgnf;
    assert(0 == (u_res & (UINT64_C(1) << 63)));
    return ((uint64_t)!!sign << 63) | u_res;
}

static inline bool
isNaNF64UI(uint64_t a)
{
    return 2047 == expF64UI(a) && 0 != fracF64UI(a);
}

static inline bool
isInf64UI(uint64_t a)
{
    return 2047 == expF64UI(a) && 0 == fracF64UI(a);
}

static inline bool
isZero64UI(uint64_t a)
{
    return 0 == (a << 1);
}

#define signExtF80UI64( a64 ) ((bool) ((uint16_t) (a64)>>15))
#define expExtF80UI64( a64 ) ((a64) & 0x7FFF)
static __inline uint16_t
packToExtF80UI64(bool sign, int32_t exp) 
{
    return (uint16_t)!!sign << 15 | (uint16_t)exp;
}

#define isNaNExtF80UI( a64, a0 ) ((((a64) & 0x7FFF) == 0x7FFF) && ((a0) & UINT64_C( 0x7FFFFFFFFFFFFFFF )))

#ifdef SOFTFLOAT_FAST_INT64
/** @bug union of same type */
union extF80M_extF80
{
    extFloat80M fM;
    extFloat80_t f;
};
/** @deprecated */
union ui128_f128
{
    uint128 ui;
    float128_t f;
};
#endif

enum
{
    softfloat_mulAdd_subC = 1,
    softfloat_mulAdd_subProd = 2
};

#ifdef SOFTFLOAT_FAST_INT64
static inline uint128
f_as_u_128(float128_t v)
{
    return *(uint128 const*)&v;
}
static inline float128_t
u_as_f_128(uint128 v)
{
    return *(float128_t const*)&v;
}
#endif  /* SOFTFLOAT_FAST_INT64 */

/**
@returns true when 32-bit unsigned integer `uiA' has the bit pattern of a
32-bit floating-point NaN.
*/
static inline bool
softfloat_isNaNF32UI(uint32_t uiA)
{
    static uint32_t exp_mask = ~(~UINT32_C(0) << 8) << 23;
    static uint32_t sgnf_mask = ~(~UINT32_C(0) << 23);
    return
        (uiA & exp_mask) == exp_mask &&
        0 != (uiA & sgnf_mask);
}
/**
@returns true when 64-bit unsigned integer `uiA' has the bit pattern of a
64-bit floating-point NaN.
*/
static inline bool
softfloat_isNaNF64UI(uint64_t uiA)
{
    static uint64_t exp_mask = ~(~UINT64_C(0) << 11) << 52;
    static uint64_t sgnf_mask = ~(~UINT64_C(0) << 52);
    return
        (uiA & exp_mask) == exp_mask &&
        0 != (uiA & sgnf_mask);
}
/**
Returns true when 32-bit unsigned integer `uiA' has the bit pattern of a
32-bit floating-point signaling NaN.
*/
static inline bool
softfloat_isSigNaNF32UI(uint32_t uiA)
{
    static uint32_t exp_mask = ~(~UINT32_C(0) << 8) << 23;
    static uint32_t quite_mask = ~(~UINT32_C(0) << 1) << 22;
    static uint32_t sgnf_mask = ~(~UINT32_C(0) << 23);
    return
        (uiA & (exp_mask | quite_mask)) == exp_mask &&
        0 != (uiA & sgnf_mask);
}

uint32_t
softfloat_roundPackToUI32(bool, uint64_t, uint8_t, bool);

#ifdef SOFTFLOAT_FAST_INT64
uint64_t
softfloat_roundPackToUI64(bool,
                          uint64_t,
                          uint64_t,
                          uint8_t,
                          bool);
#else
uint64_t
softfloat_roundPackMToUI64(bool,
                           uint32_t *,
                           uint8_t,
                           bool);
#endif

int32_t
softfloat_roundPackToI32(bool, uint64_t, uint8_t, bool);

#ifdef SOFTFLOAT_FAST_INT64
int64_t
softfloat_roundPackToI64(bool,
                         uint64_t,
                         uint64_t,
                         uint8_t,
                         bool);
#else
int64_t
softfloat_roundPackMToI64(bool,
                          uint32_t *,
                          uint8_t,
                          bool);
#endif

struct exp8_sig16
{
    int8_t exp;
    uint16_t sig;
};
exp8_sig16
    softfloat_normSubnormalF16Sig(uint16_t);

float16_t
softfloat_roundPackToF16(bool,
                         int16_t,
                         uint16_t);
float16_t
softfloat_normRoundPackToF16(bool,
                             int16_t,
                             uint16_t);

float16_t
softfloat_addMagsF16(uint16_t,
                     uint16_t);
float16_t
softfloat_subMagsF16(uint16_t,
                     uint16_t);
float16_t
softfloat_mulAddF16(uint16_t,
                    uint16_t,
                    uint16_t,
                    uint8_t);

struct exp16_sig32
{
    int16_t exp;
    uint32_t sig;
};
exp16_sig32
    softfloat_normSubnormalF32Sig(uint32_t);

float32_t
softfloat_roundPackToF32(bool,
                         int16_t,
                         uint32_t);
float32_t
softfloat_normRoundPackToF32(bool,
                             int16_t,
                             uint32_t);

float32_t
softfloat_addMagsF32(uint32_t,
                     uint32_t);
float32_t
softfloat_subMagsF32(uint32_t,
                     uint32_t);
float32_t
softfloat_mulAddF32(uint32_t,
                    uint32_t,
                    uint32_t,
                    uint8_t);

struct exp16_sig64
{
    int16_t exp;
    uint64_t sig;
};
exp16_sig64
    softfloat_normSubnormalF64Sig(uint64_t);

float64_t
softfloat_roundPackToF64(bool,
                         int16_t,
                         uint64_t);
float64_t
softfloat_normRoundPackToF64(bool,
                             int16_t,
                             uint64_t);

float64_t
softfloat_addMagsF64(uint64_t,
                     uint64_t,
                     bool);
float64_t
softfloat_subMagsF64(uint64_t,
                     uint64_t,
                     bool);
float64_t
softfloat_mulAddF64(uint64_t,
                    uint64_t,
                    uint64_t,
                    uint8_t);

#ifdef SOFTFLOAT_FAST_INT64

#define signF128UI64( a64 ) ((bool) ((uint64_t) (a64)>>63))
#define expF128UI64( a64 ) ((int32_t) ((a64)>>48) & 0x7FFF)
#define fracF128UI64( a64 ) ((a64) & UINT64_C( 0x0000FFFFFFFFFFFF ))
#define packToF128UI64( sign, exp, sig64 ) (((uint64_t) (sign)<<63) + ((uint64_t) (exp)<<48) + (sig64))

#define isNaNF128UI( a64, a0 ) (((~(a64) & UINT64_C( 0x7FFF000000000000 )) == 0) && (a0 || ((a64) & UINT64_C( 0x0000FFFFFFFFFFFF ))))

struct exp32_sig64
{
    int32_t exp;
    uint64_t sig;
};
exp32_sig64 softfloat_normSubnormalExtF80Sig(uint64_t);

extFloat80_t
softfloat_roundPackToExtF80(bool,
                            int32_t,
                            uint64_t,
                            uint64_t,
                            uint8_t);
extFloat80_t
softfloat_normRoundPackToExtF80(bool,
                                int32_t,
                                uint64_t,
                                uint64_t,
                                uint8_t);

extFloat80_t
softfloat_addMagsExtF80(uint16_t,
                        uint64_t,
                        uint16_t,
                        uint64_t,
                        bool);
extFloat80_t
softfloat_subMagsExtF80(uint16_t,
                        uint64_t,
                        uint16_t,
                        uint64_t,
                        bool);

struct exp32_sig128
{
    int32_t exp;
    uint128 sig;
};

exp32_sig128
    softfloat_normSubnormalF128Sig(uint64_t,
                                   uint64_t);

float128_t
softfloat_roundPackToF128(bool,
                          int32_t,
                          uint64_t,
                          uint64_t,
                          uint64_t);
float128_t
softfloat_normRoundPackToF128(bool,
                              int32_t,
                              uint64_t,
                              uint64_t);

float128_t
softfloat_addMagsF128(uint64_t, uint64_t, uint64_t, uint64_t, bool);
float128_t
softfloat_subMagsF128(uint64_t,
                      uint64_t,
                      uint64_t,
                      uint64_t,
                      bool);
float128_t
softfloat_mulAddF128(uint64_t,
                     uint64_t,
                     uint64_t,
                     uint64_t,
                     uint64_t,
                     uint64_t,
                     uint8_t);

#else

#define signF128UI96( a96 ) ((bool) ((uint32_t) (a96)>>31))
#define expF128UI96( a96 ) ((int32_t) ((a96)>>16) & 0x7FFF)
#define fracF128UI96( a96 ) ((a96) & 0x0000FFFF)
#define packToF128UI96( sign, exp, sig96 ) (((uint32_t) (sign)<<31) + ((uint32_t) (exp)<<16) + (sig96))

bool
softfloat_tryPropagateNaNExtF80M(extFloat80M const*,
                                 extFloat80M const*,
                                 extFloat80M*);
void
softfloat_invalidExtF80M(extFloat80M *);

int
softfloat_normExtF80SigM(uint64_t *);

void
softfloat_roundPackMToExtF80M(bool, int32_t, uint32_t *, uint8_t, extFloat80M *);
void
softfloat_normRoundPackMToExtF80M(bool, int32_t, uint32_t *, uint8_t, extFloat80M *);

void
softfloat_addExtF80M(extFloat80M const*,
                     extFloat80M const*,
                     extFloat80M*,
                     bool);

int
softfloat_compareNonnormExtF80M(extFloat80M const*,
                                extFloat80M const*);

bool softfloat_isNaNF128M(uint32_t const*);

bool
softfloat_tryPropagateNaNF128M(uint32_t const*,
                               uint32_t const*,
                               uint32_t *);

void
softfloat_invalidF128M(uint32_t *);

int
softfloat_shiftNormSigF128M(uint32_t const *, 
                            uint8_t, 
                            uint32_t *);

void
softfloat_roundPackMToF128M(bool, 
                            int32_t, 
                            uint32_t *, 
                            uint32_t *);
void
softfloat_normRoundPackMToF128M(bool, 
                                int32_t, 
                                uint32_t *, 
                                uint32_t *);

void
softfloat_addF128M(uint32_t const*,
                   uint32_t const *, 
                   uint32_t *, 
                   bool);
void
softfloat_mulAddF128M(uint32_t const*,
                      uint32_t const*,
                      uint32_t const*,
                      uint32_t*,
                      uint8_t);

#endif

#endif  /* SOFTFLOAT_INTERNALS_H_ */
