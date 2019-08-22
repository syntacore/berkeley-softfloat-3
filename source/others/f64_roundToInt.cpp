
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

float64_t
f64_roundToInt(float64_t const a,
               uint8_t const roundingMode,
               bool const exact)
{
    using namespace softfloat::internals;
    uint64_t const uiA = f_as_u(a);
    int16_t const exp = get_exp(uiA);

    if (exp <= 0x3FE) {
        if (0 == (uiA & INT64_MAX)) {
            return a;
        }

        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        uint64_t const uiZ = uiA & packToF64UI(1, 0, 0);

        switch (roundingMode) {
        case softfloat_round_near_even:
            if (!get_frac(uiA)) {
                return u_as_f(uiZ);
                break;
            }

            [[fallthrough]];

        case softfloat_round_near_maxMag:
            if (exp == 0x3FE) {
                return u_as_f(uiZ | packToF64UI(0, 0x3FF, 0));
            }

            break;

        case softfloat_round_min:
            if (uiZ) {
                return u_as_f(packToF64UI(1, 0x3FF, 0));
            }

            break;

        case softfloat_round_max:
            if (!uiZ) {
                return u_as_f(packToF64UI(0, 0x3FF, 0));
            }

            break;
        }

        return u_as_f(uiZ);
    }

    if (0x433 <= exp) {
        if (exp == 0x7FF && get_frac(uiA)) {
            return u_as_f(propagate_NaN(uiA, 0));
        }

        return a;
    }

    uint64_t uiZ = uiA;
    uint64_t const lastBitMask = static_cast<uint64_t>(1) << (0x433 - exp);
    uint64_t const roundBitsMask = lastBitMask - 1;

    if (softfloat_round_near_maxMag == roundingMode) {
        uiZ += lastBitMask >> 1;
    } else if (softfloat_round_near_even == roundingMode) {
        uiZ += lastBitMask >> 1;

        if (0 == (uiZ & roundBitsMask)) {
            uiZ &= ~lastBitMask;
        }
    } else if (softfloat_round_minMag != roundingMode && is_sign(uiZ) != (softfloat_round_max == roundingMode)) {
        uiZ += roundBitsMask;
    }

    uiZ &= ~roundBitsMask;

    if (exact && uiZ != uiA) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return u_as_f(uiZ);
}
