
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

float128_t
f16_to_f128(float16_t const a)
{
    using namespace softfloat::internals;
    uint16_t const uiA = f_as_u_16(a);
    bool const sign = is_sign(uiA);
    int8_t const exp = expF16UI(uiA);
    uint16_t const frac = fracF16UI(uiA);

    if (0x1F == exp) {
        // infinity or NaN
        return
            float128_t(
                0 != frac ?
                softfloat_commonNaNToF128UI(softfloat_f16UIToCommonNaN(uiA)) :
                uint128{packToF128UI64(sign, 0x7FFF, 0), 0});
    }

    assert(0x1F != exp);

    if (0 != exp) {
        // normalized
        return
            float128_t(
                uint128(
                    packToF128UI64(sign,
                                   exp + 0x3FF0,
                                   static_cast<uint64_t>(frac) << 38),
                    0));
    }

    assert(0 == exp);

    if (0 == frac) {
        assert(0 == exp && 0 == frac);
        return float128_t(uint128{packToF128UI64(sign, 0, 0), 0});
    }

    assert(0 == exp && 0 != frac);

    exp8_sig16 const normExpSig(frac);
    return
        float128_t(
            uint128(
                packToF128UI64(sign,
                               normExpSig.exp - 1 + 0x3FF0,
                               static_cast<uint64_t>(normExpSig.sig) << 38),
                0));
}
