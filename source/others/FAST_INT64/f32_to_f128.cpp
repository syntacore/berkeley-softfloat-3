
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

#include "softfloat/functions.h"

#include "internals.hpp"
#include "specialize.hpp"

float128_t
f32_to_f128(float32_t a)
{
    using namespace softfloat;
    uint32_t const uiA = f_as_u_32(a);
    bool const sign = signF32UI(uiA);
    int16_t exp = expF32UI(uiA);
    uint32_t frac = fracF32UI(uiA);

    if (exp == 0xFF) {
        if (frac) {
            return u_as_f_128(softfloat_commonNaNToF128UI(softfloat_f32UIToCommonNaN(uiA)));
        }

        uint128 uiZ;
        uiZ.v64 = packToF128UI64(sign, 0x7FFF, 0);
        uiZ.v0 = 0;
        return u_as_f_128(uiZ);
    }

    if (0 == exp) {
        if (0 == frac) {
            struct uint128 uiZ;
            uiZ.v64 = packToF128UI64(sign, 0, 0);
            uiZ.v0 = 0;
            return u_as_f_128(uiZ);
        }

        exp16_sig32 const normExpSig = softfloat_normSubnormalF32Sig(frac);
        exp = normExpSig.exp - 1;
        frac = normExpSig.sig;
    }

    uint128 uiZ;
    uiZ.v64 = packToF128UI64(sign, exp + 0x3F80, (uint64_t)frac << 25);
    uiZ.v0 = 0;
    return u_as_f_128(uiZ);
}

