
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
All rights reserved.

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

#ifdef _MSC_VER
#pragma warning( disable : 5030)
#endif

float32_t
f32_roundToInt(float32_t const a,
               uint8_t const roundingMode,
               bool const exact)
{
    using namespace softfloat::internals;
    uint32_t const uiA = f_as_u(a);
    int16_t const exp = get_exp(a);

    if (exp <= 0x7E) {
        if (0 == static_cast<uint32_t>(uiA << 1)) {
            // is zero
            return a;
        }

        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        uint32_t const uiZ = uiA & make_signed_zero<uint32_t>(true);

        switch (roundingMode) {
        case softfloat_round_near_even:
            if (0 == get_frac(a)) {
                break;
            }

            [[fallthrough]];

        case softfloat_round_near_maxMag:
            if (0x7E == exp) {
                return u_as_f(uiZ | pack_to_F32_UI(false, 0x7F, 0));
            }

            break;

        case softfloat_round_min:
            if (0 != uiZ) {
                return u_as_f(pack_to_F32_UI(true, 0x7F, 0));
            }

            break;

        case softfloat_round_max:
            if (0 == uiZ) {
                return u_as_f(pack_to_F32_UI(false, 0x7F, 0));
            }

            break;
        }

        return u_as_f(uiZ);
    }

    if (0x96 <= exp) {
        return is_NaN(a) ? propagate_NaN(a, a) : a;
    }

    uint32_t uiZ = uiA;
    uint32_t const lastBitMask = static_cast<uint32_t>(1) << (0x96 - exp);
    uint32_t const roundBitsMask = lastBitMask - 1;

    if (softfloat_round_near_maxMag == roundingMode) {
        uiZ += lastBitMask >> 1;
    } else if (roundingMode == softfloat_round_near_even) {
        uiZ += lastBitMask >> 1;

        if (0 == (uiZ & roundBitsMask)) {
            uiZ &= ~lastBitMask;
        }
    } else if (softfloat_round_minMag != roundingMode && (is_sign(uiZ) != (softfloat_round_max == roundingMode))) {
        uiZ += roundBitsMask;
    }

    uiZ &= ~roundBitsMask;

    if (exact && uiZ != uiA) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return u_as_f(uiZ);
}
