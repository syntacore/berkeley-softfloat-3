
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

void
f32_to_f128M(float32_t a,
             float128_t* zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    *zPtr = f32_to_f128(a);
#else
    using namespace softfloat::internals::slow_int64;
    uint32_t* zWPtr = reinterpret_cast<uint32_t*>(zPtr);
    bool const sign = is_sign(a);
    int16_t exp = get_exp(a);
    uint32_t frac = get_frac(a);

    if (is_NaN(a)) {
        commonNaN_to_F128M(commonNaN_from_f32UI(f_as_u(a)), zWPtr);
        return;
    }

    if (is_inf(a)) {
        zWPtr[indexWord(4, 3)] = packToF128UI96(sign, 0x7FFF, 0);
        zWPtr[indexWord(4, 2)] = 0;
        zWPtr[indexWord(4, 1)] = 0;
        zWPtr[indexWord(4, 0)] = 0;
        return;
    }

    if (is_zero(a)) {
        zWPtr[indexWord(4, 3)] = packToF128UI96(sign, 0, 0);
        zWPtr[indexWord(4, 2)] = 0;
        zWPtr[indexWord(4, 1)] = 0;
        zWPtr[indexWord(4, 0)] = 0;
        return;
    }

    if (0 == exp) {
        exp16_sig32 const normExpSig(frac);
        exp = normExpSig.exp - 1;
        frac = normExpSig.sig;
    }

    uint64_t const frac64 = static_cast<uint64_t>(frac) << 25;
    zWPtr[indexWord(4, 3)] = packToF128UI96(sign, exp + 0x3F80u, frac64 >> 32);
    zWPtr[indexWord(4, 2)] = static_cast<uint32_t>(frac64);
    zWPtr[indexWord(4, 1)] = 0;
    zWPtr[indexWord(4, 0)] = 0;
    return;
#endif
}
