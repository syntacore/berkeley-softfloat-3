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

namespace softfloat {
namespace internals {

inline constexpr bool
signF16UI(uint16_t a)
{
    return static_cast<int16_t>(a) < 0;
}

inline constexpr int8_t
expF16UI(uint16_t a)
{
    return static_cast<int8_t>((a >> 10) & 0x1F);
}

inline constexpr uint16_t
fracF16UI(uint16_t a)
{
    return a & 0x03FFu;
}

inline constexpr uint16_t
packToF16UI(bool sign,
            int8_t exp,
            uint16_t sig)
{
    return
        static_cast<uint16_t>(
        (static_cast<uint16_t>(!!sign) << 15) +
        (static_cast<uint16_t>(exp) << 10) +
            sig);
}

inline constexpr bool
isNaNF16UI(uint16_t a)
{
    return
        0 == (UINT16_C(0x7C00) & ~a) &&
        0 != (UINT16_C(0x03FF) &  a);
}

inline constexpr uint16_t
f_as_u_16(float16_t v)
{
    return v.v;
}

inline constexpr float16_t
u_as_f_16(uint16_t v)
{
    return float16_t{v};
}

inline constexpr uint32_t
f_as_u_32(float32_t v)
{
    return v.v;
}

inline constexpr float32_t
u_as_f_32(uint32_t v)
{
    return float32_t{v};
}

inline constexpr uint64_t
f_as_u_64(float64_t v)
{
    return v.v;
}

inline constexpr float64_t
u_as_f_64(uint64_t v)
{
    return float64_t{v};
}

inline constexpr bool
signF32UI(uint32_t a)
{
    return static_cast<int32_t>(a) < 0;
}

/** @bug return signed 16-bits value instead unsigned 8-bits value */
inline constexpr int16_t
expF32UI(uint32_t a)
{
    return static_cast<int16_t>((a >> 23) & ~(~UINT32_C(0) << 8));
}

inline constexpr uint32_t
fracF32UI(uint32_t a)
{
    return a & ~(~UINT32_C(0) << 23);
}

inline constexpr uint32_t
packToF32UI(bool sign,
            int16_t exp,
            uint32_t sgnf)
{
    return
        (static_cast<uint32_t>(!!sign) << 31) |
        ((static_cast<uint32_t>(exp) << 23) + sgnf);
}

inline constexpr uint32_t
signed_zero_F32UI(bool sign)
{
    return packToF32UI(sign, 0, 0u);
}

inline constexpr float32_t
signed_zero_F32(bool sign)
{
    return u_as_f_32(signed_zero_F32UI(sign));
}

inline constexpr uint32_t
signed_inf_F32UI(bool sign)
{
    return packToF32UI(sign, 0xFF, 0u);
}

inline constexpr float32_t
signed_inf_F32(bool sign)
{
    return u_as_f_32(signed_inf_F32UI(sign));
}

inline constexpr bool
isNaNF32UI(uint32_t a)
{
    return
        255 == expF32UI(a) &&
        0 != fracF32UI(a);
}

inline constexpr bool
isInf32UI(uint32_t a)
{
    return
        255 == expF32UI(a) &&
        0 == fracF32UI(a);
}

inline constexpr bool
isZero32UI(uint32_t a)
{
    return
        0 == expF32UI(a) &&
        0 == fracF32UI(a);
}

inline constexpr bool
isSubnormal32UI(uint32_t a)
{
    return
        0 == expF32UI(a) &&
        0 != fracF32UI(a);
}

inline constexpr bool
signF64UI(uint64_t a)
{
    return static_cast<int64_t>(a) < 0;
}

inline constexpr int16_t
expF64UI(uint64_t a)
{
    return static_cast<int16_t>((a >> 52) & ~(~UINT32_C(0) << 11));
}

inline constexpr uint64_t
fracF64UI(uint64_t a)
{
    return a & ~(~UINT64_C(0) << 52);
}

inline constexpr uint64_t
packToF64UI(bool sign, int16_t exp, uint64_t sgnf)
{
    return
        (static_cast<uint64_t>(!!sign) << 63) |
        ((static_cast<uint64_t>(exp) << 52) + sgnf);
}

inline constexpr bool
isNaNF64UI(uint64_t a)
{
    return
        2047 == expF64UI(a) &&
        0 != fracF64UI(a);
}

inline constexpr bool
isInf64UI(uint64_t a)
{
    return (~(~UINT64_C(0) << 11) << 52) == (~(~UINT64_C(0) << 63) & a);
}

inline constexpr bool
isZero64UI(uint64_t a)
{
    return 0 == (~(~UINT64_C(0) << 63) & a);
}

inline constexpr bool
signExtF80UI64(uint16_t a64)
{
    return 0 != (a64 >> 15);
}

inline constexpr uint16_t
expExtF80UI64(uint16_t a64)
{
    return static_cast<uint16_t>(0x7FFF & a64);
}

inline constexpr uint16_t
packToExtF80UI64(bool sign,
                 uint16_t exp)
{
    return static_cast<uint16_t>((!!sign << 15) | exp);
}

inline constexpr bool
isNaNExtF80UI(uint16_t a64,
              uint64_t a0)
{
    return
        UINT16_C(0x7FFF) == (UINT16_C(0x7FFF) & a64) &&
        0 != (UINT64_C(0x7FFFFFFFFFFFFFFF) & a0);
}

#ifdef SOFTFLOAT_FAST_INT64
/** @deprecated */
union ui128_f128
{
    uint128 ui;
    float128_t f;
};

inline uint128 const&
f_as_u_128(float128_t const& v)
{
    return reinterpret_cast<uint128 const&>(v);
}

inline float128_t const&
u_as_f_128(uint128 const& v)
{
    return reinterpret_cast<float128_t const&>(v);
}

#endif  /* SOFTFLOAT_FAST_INT64 */

/**
@returns true when 32-bit unsigned integer `uiA' has the bit pattern of a
32-bit floating-point NaN.
*/
inline constexpr bool
softfloat_isNaNF32UI(uint32_t uiA)
{
    return
        (~(~UINT32_C(0) << 8) << 23) == (uiA & (~(~UINT32_C(0) << 8) << 23)) &&
        0 != (uiA & (~(~UINT32_C(0) << 23)));
}

/**
@returns true when 64-bit unsigned integer `uiA' has the bit pattern of a
64-bit floating-point NaN.
*/
inline constexpr bool
softfloat_isNaNF64UI(uint64_t uiA)
{
    return
        (~(~UINT64_C(0) << 11) << 52) == ((~(~UINT64_C(0) << 11) << 52) & uiA) &&
        0 != ((~(~UINT64_C(0) << 52)) & uiA);
}
/**
Returns true when 32-bit unsigned integer `uiA' has the bit pattern of a
32-bit floating-point signaling NaN.
*/
inline constexpr bool
softfloat_isSigNaNF32UI(uint32_t uiA)
{
    return
        (~(~UINT32_C(0) << 8) << 23) == (((~(~UINT32_C(0) << 8) << 23) | (~(~UINT32_C(0) << 1) << 22)) & uiA) &&
        0 != (~(~UINT32_C(0) << 23) & uiA);
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
                           uint32_t*,
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
                          uint32_t*,
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
enum Mul_add_operations
{
    softfloat_mulAdd_madd = 0,
    softfloat_mulAdd_subC = 1,
    softfloat_mulAdd_subProd = 2
};

float16_t
softfloat_mulAddF16(uint16_t,
                    uint16_t,
                    uint16_t,
                    Mul_add_operations);

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

#ifdef SOFTFLOAT_FAST_INT64

inline constexpr bool
signF128UI64(uint64_t a64)
{
    return static_cast<int64_t>(a64) < 0;
}

inline constexpr int32_t
expF128UI64(uint64_t a64)
{
    return static_cast<int32_t>(a64 >> 48) & INT32_C(0x7FFF);
}

inline constexpr uint64_t
fracF128UI64(uint64_t a64)
{
    return a64 & UINT64_C(0x0000FFFFFFFFFFFF);
}

inline constexpr uint64_t
packToF128UI64(bool sign,
               int32_t exp,
               uint64_t sig64)
{
    return
        (static_cast<uint64_t>(sign) << 63) + 
        (static_cast<uint64_t>(exp) << 48) + 
        sig64;
}

inline constexpr bool
isNaNF128UI(uint64_t a64,
            uint64_t a0)
{
    return
        0 == (~a64 & UINT64_C(0x7FFF000000000000)) &&
        (0 != a0 || 0 != (a64 & UINT64_C(0x0000FFFFFFFFFFFF)));
}

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
                     Mul_add_operations);

#else

inline constexpr bool
signF128UI96(uint32_t a96)
{
    return static_cast<int32_t>(a96) < 0;
}

inline constexpr uint16_t
expF128UI96(uint32_t a96)
{
    return static_cast<uint16_t>((a96 >> 16) & UINT16_C(0x7FFF));
}

inline constexpr uint32_t
fracF128UI96(uint32_t a96)
{
    return a96 & 0x0000FFFF;
}

inline constexpr uint32_t
packToF128UI96(bool sign,
               unsigned exp,
               uint32_t sig96)
{
    return
        sig96 +
        (static_cast<uint32_t>((!!sign) << 31) + 
        (static_cast<uint32_t>(exp) << 16));
}

bool
softfloat_tryPropagateNaNExtF80M(extFloat80M const*,
                                 extFloat80M const*,
                                 extFloat80M*);
void
softfloat_invalidExtF80M(extFloat80M*);

int
softfloat_normExtF80SigM(uint64_t*);

void
softfloat_roundPackMToExtF80M(bool, int32_t, uint32_t*, uint8_t, extFloat80M*);
void
softfloat_normRoundPackMToExtF80M(bool, int32_t, uint32_t*, uint8_t, extFloat80M*);

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
                               uint32_t*);

void
softfloat_invalidF128M(uint32_t*);

int
softfloat_shiftNormSigF128M(uint32_t const*,
                            uint8_t,
                            uint32_t*);

void
softfloat_roundPackMToF128M(bool,
                            int32_t,
                            uint32_t*,
                            uint32_t*);
void
softfloat_normRoundPackMToF128M(bool,
                                int32_t,
                                uint32_t*,
                                uint32_t*);

void
softfloat_addF128M(uint32_t const*,
                   uint32_t const*,
                   uint32_t*,
                   bool);
void
softfloat_mulAddF128M(uint32_t const*,
                      uint32_t const*,
                      uint32_t const*,
                      uint32_t*,
                      Mul_add_operations);

#endif
}  // namespace internals
}  // namespace softfloat

#endif  /* SOFTFLOAT_INTERNALS_H_ */
