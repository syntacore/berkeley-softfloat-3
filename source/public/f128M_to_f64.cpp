
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

float64_t
f128M_to_f64(const float128_t* aPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    return f128_to_f64(*aPtr);
#else
    using namespace softfloat::internals::slow_int64;

    uint32_t const* aWPtr = reinterpret_cast<const uint32_t*>(aPtr);
    uint32_t const uiA96 = aWPtr[indexWordHi(4)];
    bool const sign = is_sign(uiA96);
    int32_t exp = expF128UI96(uiA96);
    uint64_t frac64 = static_cast<uint64_t>(fracF128UI96(uiA96)) << 32 | aWPtr[indexWord(4, 2)];

    if (exp == 0x7FFF) {
        if (frac64 || aWPtr[indexWord(4, 1)] | aWPtr[indexWord(4, 0)]) {
            return u_as_f(commonNaN_to_F64UI(commonNaN_from_f128M(aWPtr)));
        }

        return u_as_f(packToF64UI(sign, 0x7FF, 0));
    }

    uint32_t const frac32 = aWPtr[indexWord(4, 1)];
    frac64 = frac64 << 14 | frac32 >> 18;

    if (0 != (frac32 & 0x0003FFFF) || 0 != aWPtr[indexWord(4, 0)]) {
        frac64 |= 1;
    }

    if (0 == exp && 0 == frac64) {
        return u_as_f(packToF64UI(sign, 0, 0));
    }

    exp -= 0x3C01;

    if (exp < -0x1000) {
        exp = -0x1000;
    }

    /**
    @todo Warning   C4242   'function': conversion from 'int32_t' to 'int16_t', possible loss of data
    */
    return
        roundPackToF64(sign, static_cast<int16_t>(exp), frac64 | UINT64_C(0x4000000000000000));
#endif
}
