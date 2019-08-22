
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015 The Regents of the University of
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

#include "model.hpp"

#ifndef SOFTFLOAT_FAST_INT64
#error Fast int64_t operations only
#endif

namespace softfloat {
namespace internals {
namespace {

/**
Returns true if the 128-bit unsigned integer formed by concatenating `a64'
and `a0' is equal to the 128-bit unsigned integer formed by concatenating
`b64' and `b0'.
*/
static inline constexpr bool
softfloat_eq128(uint64_t const& a64,
                uint64_t const& a0,
                uint64_t const& b64,
                uint64_t const& b0)
{
    return a64 == b64 && a0 == b0;
}

static inline constexpr bool
softfloat_eq128(uint128 const& a,
                uint128 const& b)
{
    return softfloat_eq128(a.v64, a.v0, b.v64, b.v0);
}

}  // namespace

float128_t
softfloat_roundPackToF128(bool const sign,
                          int32_t exp,
                          uint64_t sig64,
                          uint64_t sig0,
                          uint64_t sigExtra)
{
    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
    bool const roundNearEven = softfloat_round_near_even == softfloat_roundingMode;
    bool doIncrement = UINT64_C(0x8000000000000000) <= sigExtra;

    if (!roundNearEven && softfloat_round_near_maxMag != softfloat_roundingMode) {
        doIncrement =
            (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode &&
            0 != sigExtra;
    }

    if (0x7FFD <= static_cast<uint32_t>(exp)) {
        if (exp < 0) {
            bool const isTiny =
                softfloat_tininess_beforeRounding == softfloat_detectTininess ||
                exp < -1 ||
                !doIncrement ||
                softfloat_lt128(uint128(sig64, sig0), uint128(UINT64_C(0x0001FFFFFFFFFFFF), UINT64_MAX));

            uint128_extra const sig128Extra = softfloat_shiftRightJam128Extra(uint128_extra(uint128(sig64, sig0), sigExtra), static_cast<uint32_t>(-exp));
            sig64 = sig128Extra.v.v64;
            sig0 = sig128Extra.v.v0;
            sigExtra = sig128Extra.extra;
            exp = 0;

            if (isTiny && 0 != sigExtra) {
                softfloat_raiseFlags(softfloat_flag_underflow);
            }

            doIncrement = UINT64_C(0x8000000000000000) <= sigExtra;

            if (!roundNearEven && softfloat_round_near_maxMag != softfloat_roundingMode) {
                doIncrement =
                    (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode &&
                    0 != sigExtra;
            }
        } else if (
            0x7FFD < exp ||
            (
                0x7FFD == exp &&
                doIncrement &&
                softfloat_eq128(uint128(sig64, sig0), uint128(UINT64_C(0x0001FFFFFFFFFFFF), UINT64_C(0xFFFFFFFFFFFFFFFF)))
            )
        ) {
            softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);

            if (
                roundNearEven ||
                softfloat_round_near_maxMag == softfloat_roundingMode ||
                (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode
            ) {
                return float128_t(uint128{packToF128UI64(sign, 0x7FFF, 0), 0});
            }

            return float128_t(uint128{packToF128UI64(sign, 0x7FFE, UINT64_C(0x0000FFFFFFFFFFFF)), UINT64_MAX});
        }
    }

    if (0 != sigExtra) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    if (doIncrement) {
        uint128 const sig128 = softfloat_add128(uint128{sig64, sig0}, uint128{0, 1});
        return float128_t(uint128{
            packToF128UI64(sign, exp, sig128.v64),
            sig128.v0 & ~static_cast<uint64_t>(!(sigExtra & INT64_MAX) & roundNearEven)});
    }

    if (0 == (sig64 | sig0)) {
        exp = 0;
    }

    return float128_t(uint128{packToF128UI64(sign, exp, sig64), sig0});
}

}  // namespace internals
}  // namespace softfloat
