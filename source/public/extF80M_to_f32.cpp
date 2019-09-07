
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015 The Regents of the University of
California.  All rights reserved.
*/
/*
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

float32_t
extF80M_to_f32(const extFloat80_t* aPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    return extF80_to_f32(*aPtr);
#else
    using namespace softfloat::internals::slow_int64;

    extFloat80M const* aSPtr = aPtr;
    uint16_t const uiA64 = aSPtr->signExp;
    bool const sign = is_sign(uiA64);
    int32_t exp = expExtF80UI64(uiA64);
    uint64_t sig = aSPtr->signif;

    if (0x7FFF == exp) {
        if (0 != (sig & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
            return u_as_f(softfloat_commonNaNToF32UI(commonNaN{*aSPtr}));
        }

        return make_signed_inf<float32_t>(sign);
    }

    if (0 == (sig & UINT64_C(0x8000000000000000))) {
        if (0 == sig) {
            return make_signed_zero<float32_t>(sign);
        }

        exp += softfloat_normExtF80SigM(&sig);
    }

    uint32_t const sig32 = static_cast<uint32_t>(softfloat_shortShiftRightJam64(sig, 33));
    exp -= 0x3F81;

    if (exp < -0x1000) {
        exp = -0x1000;
    }

    assert(INT16_MIN <= exp && exp <= INT16_MAX);
    return softfloat_roundPackToF32(sign, static_cast<int16_t>(exp), sig32);
#endif
}
