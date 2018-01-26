
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

#include "target.hpp"

float16_t
f16_roundToInt(float16_t a,
               uint8_t roundingMode,
               bool exact)
{
    using namespace softfloat::internals;
    uint16_t const uiA = f_as_u_16(a);
    int8_t const exp = expF16UI(uiA);

    if (exp <= 0xE) {
        if (!static_cast<uint16_t>(uiA << 1)) {
            return a;
        }

        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        uint16_t const uiZ = static_cast<uint16_t>(uiA & packToF16UI(1, 0, 0u));

        switch (roundingMode) {
        case softfloat_round_near_even:
            if (!fracF16UI(uiA)) {
                return u_as_f_16(uiZ);
            }

        case softfloat_round_near_maxMag:
            if (exp == 0xE) {
                return u_as_f_16(static_cast<uint16_t>(uiZ | packToF16UI(0, 0xF, 0u)));
            }

            break;

        case softfloat_round_min:
            if (uiZ) {
                return u_as_f_16(packToF16UI(1, 0xF, 0));
            }

            break;

        case softfloat_round_max:
            if (!uiZ) {
                return u_as_f_16(packToF16UI(0, 0xF, 0));
            }

            break;
        }

        return u_as_f_16(uiZ);
    } else if (0x19 <= exp) {
        return
            exp == 0x1F && 0 != fracF16UI(uiA) ?
            u_as_f_16(softfloat_propagateNaNF16UI(uiA, 0)) : a;
    } else {
        uint16_t uiZ = uiA;
        uint16_t const lastBitMask = static_cast<uint16_t>(1 << (0x19 - exp));
        uint16_t const roundBitsMask = lastBitMask - 1u;

        if (roundingMode == softfloat_round_near_maxMag) {
            uiZ += lastBitMask >> 1;
        } else if (roundingMode == softfloat_round_near_even) {
            uiZ += lastBitMask >> 1;

            if (!(uiZ & roundBitsMask)) {
                uiZ &= ~lastBitMask;
            }
        } else if (roundingMode != softfloat_round_minMag) {
            if (signF16UI(uiZ) ^ (roundingMode == softfloat_round_max)) {
                uiZ += roundBitsMask;
            }
        }

        uiZ &= ~roundBitsMask;

        if (exact && (uiZ != uiA)) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return u_as_f_16(uiZ);
    }
}
