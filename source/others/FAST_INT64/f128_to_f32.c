
/*============================================================================

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

=============================================================================*/

#include "softfloat/functions.h"

#include "internals.h"
#include "specialize.h"

#include <assert.h>

float32_t f128_to_f32(float128_t a)
{
    union ui128_f128 uA;
    uint64_t uiA64, uiA0;
    bool sign;
    int32_t exp;
    uint64_t frac64;
    uint32_t uiZ;
    union ui32_f32 uZ;

    uA.f = a;
    uiA64 = uA.ui.v64;
    uiA0 = uA.ui.v0;
    sign = signF128UI64(uiA64);
    exp = expF128UI64(uiA64);
    frac64 = fracF128UI64(uiA64) | (uiA0 != 0);
    if (exp != INT16_MAX) {
        uint32_t const frac32 = softfloat_shortShiftRightJam64(frac64, 18);
        if (exp | frac32) {
            exp -= 0x3F81;
            if (sizeof(int16_t) < sizeof exp) {
                if (exp < -0x1000) {exp = -0x1000;}
            }
            assert(INT16_MIN <= exp && exp <= INT16_MAX);
            return softfloat_roundPackToF32(sign, (int16_t)exp, frac32 | 0x40000000);
        } else {
            uiZ = packToF32UI(sign, 0, 0);
        }
    } else {
        if (frac64) {
            struct commonNaN commonNaN;
            softfloat_f128UIToCommonNaN(uiA64, uiA0, &commonNaN);
            uiZ = softfloat_commonNaNToF32UI(&commonNaN);
        } else {
            uiZ = packToF32UI(sign, 0xFF, 0);
        }
    }
    uZ.ui = uiZ;
    return uZ.f;

}

