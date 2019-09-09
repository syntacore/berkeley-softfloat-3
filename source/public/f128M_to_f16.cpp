
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

float16_t
f128M_to_f16(const float128_t* aPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    return f128_to_f16(*aPtr);
#else
    using namespace softfloat::internals::slow_int64;
    uint32_t const* const aWPtr = reinterpret_cast<uint32_t const*>(aPtr);

    uint32_t const uiA96 = aWPtr[index_word_hi(4)];
    bool const sign = is_sign(uiA96);
    int32_t exp = exp_F128_UI96(uiA96);
    uint32_t const frac32 =
        frac_F128_UI96(uiA96) |
        ((aWPtr[index_word(4, 2)] | aWPtr[index_word(4, 1)] | aWPtr[index_word(4, 0)]) != 0);

    if (exp == 0x7FFF) {
        if (frac32) {
            return u_as_f(commonNaN_to_F16UI(commonNaN_from_M_f128(aWPtr)));
        }

        return u_as_f(pack_to_F16_UI(sign, 0x1F, 0));
    }

    uint16_t const frac16 = static_cast<uint16_t>(frac32 >> 2 | (frac32 & 3));

    if (0 == (exp | frac16)) {
        return u_as_f(pack_to_F16_UI(sign, 0, 0));
    }

    exp -= 0x3FF1;

    if (exp < -0x40) {
        exp = -0x40;
    }

    return round_pack_to_F16(sign, static_cast<int16_t>(exp), static_cast<uint16_t>(frac16 | 0x4000));
#endif
}
