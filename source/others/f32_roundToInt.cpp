
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

#include "softfloat/functions.h"

#include "internals.hpp"
#include "specialize.hpp"

float32_t
f32_roundToInt(float32_t a,
               uint8_t roundingMode,
               bool exact)
{
    using namespace softfloat::internals;
    uint32_t const uiA = f_as_u_32(a);
    int16_t const exp = expF32UI(uiA);

    if (exp <= 0x7E) {
        if (!static_cast<uint32_t>(uiA << 1)) {
            return a;
        }

        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        uint32_t const uiZ = uiA & signed_zero_F32UI(true);

        switch (roundingMode) {
        case softfloat_round_near_even:
            if (!fracF32UI(uiA)) {
                break;
            }

        case softfloat_round_near_maxMag:
            if (exp == 0x7E) {
                return u_as_f_32(uiZ | packToF32UI(false, 0x7F, 0));
            }

            break;

        case softfloat_round_min:
            if (uiZ) {
                return u_as_f_32(packToF32UI(true, 0x7F, 0));
            }

            break;

        case softfloat_round_max:
            if (!uiZ) {
                return u_as_f_32(packToF32UI(false, 0x7F, 0));
            }

            break;
        }

        return u_as_f_32(uiZ);
    }

    if (0x96 <= exp) {
        return
            exp == 0xFF && fracF32UI(uiA) ? u_as_f_32(softfloat_propagateNaNF32UI(uiA, 0)) : a;
    }

    uint32_t uiZ = uiA;
    uint32_t const lastBitMask = static_cast<uint32_t>(1) << (0x96 - exp);
    uint32_t const roundBitsMask = lastBitMask - 1;

    if (roundingMode == softfloat_round_near_maxMag) {
        uiZ += lastBitMask >> 1;
    } else if (roundingMode == softfloat_round_near_even) {
        uiZ += lastBitMask >> 1;

        if (!(uiZ & roundBitsMask)) {
            uiZ &= ~lastBitMask;
        }
    } else if (roundingMode != softfloat_round_minMag) {
        if (signF32UI(uiZ) ^ (roundingMode == softfloat_round_max)) {
            uiZ += roundBitsMask;
        }
    }

    uiZ &= ~roundBitsMask;

    if (exact && (uiZ != uiA)) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return u_as_f_32(uiZ);
}
