
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
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

#if (SOFTFLOAT_FAST_INT64)

int64_t
f128M_to_i64_r_minMag(const float128_t* aPtr,
                      bool exact)
{
    return f128_to_i64_r_minMag(*aPtr, exact);
}

#else

namespace {

static inline int64_t
invalid(bool const sign,
        int32_t const exp,
        uint32_t sig96,
        uint32_t const* aWPtr)
{
    using namespace softfloat::internals;
    softfloat_raiseFlags(softfloat_flag_invalid);
    return
        exp == 0x7FFF && (0 != sig96 || 0 != (aWPtr[index_word(4, 2)] | aWPtr[index_word(4, 1)] | aWPtr[index_word(4, 0)])) ?
        i64_fromNaN :
        sign ?
        i64_fromNegOverflow :
        i64_fromPosOverflow;
}

}  // namespace

int64_t
f128M_to_i64_r_minMag(const float128_t* aPtr,
                      bool exact)
{
    using namespace softfloat::internals::slow_int64;

    uint32_t sig[4];
    uint64_t uiZ;
    auto const aWPtr = reinterpret_cast<const uint32_t*>(aPtr);
    uint32_t const uiA96 = aWPtr[index_word_hi(4)];
    bool const sign = is_sign(uiA96);
    int32_t const exp = exp_F128_UI96(uiA96);
    uint32_t sig96 = frac_F128_UI96(uiA96);

    int32_t const shiftDist = 0x403E - exp;

    if (shiftDist < 0) {
        return invalid(sign, exp, sig96, aWPtr);
    }

    if (exact) {
        if (exp) {
            sig96 |= 0x00010000;
        }

        sig[index_word(4, 3)] = sig96;
        sig[index_word(4, 2)] = aWPtr[index_word(4, 2)];
        sig[index_word(4, 1)] = aWPtr[index_word(4, 1)];
        sig[index_word(4, 0)] = aWPtr[index_word(4, 0)];
        shift_right_jam_M_128(sig, static_cast<uint8_t>(shiftDist + 17), sig);
        uiZ = static_cast<uint64_t>(sig[index_word(4, 2)]) << 32 | sig[index_word(4, 1)];

        if (uiZ >> 63 && (!sign || (uiZ != UINT64_C(0x8000000000000000)))) {
            return invalid(sign, exp, sig96, aWPtr);
        }

        if (sig[index_word_lo(4)]) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }
    } else {
        if (64 <= shiftDist) {
            return 0;
        }

        uiZ =
            static_cast<uint64_t>(sig96) << 47 |
            static_cast<uint64_t>(aWPtr[index_word(4, 2)]) << 15 |
            aWPtr[index_word(4, 1)] >> 17;

        if (shiftDist) {
            uiZ |= UINT64_C(0x8000000000000000);
            uiZ >>= shiftDist;
        } else {
            if (uiZ || !sign) {
                return invalid(sign, exp, sig96, aWPtr);
            }

            uiZ |= UINT64_C(0x8000000000000000);
        }
    }

    return sign ? -static_cast<int64_t>(uiZ) : static_cast<int64_t>(uiZ);
}
#endif
