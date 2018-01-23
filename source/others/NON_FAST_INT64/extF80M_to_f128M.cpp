
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
All rights reserved.
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

#include "softfloat/functions.h"

#include "internals.hpp"
#include "target.hpp"

void
extF80M_to_f128M(extFloat80_t const* const aPtr, float128_t* const zPtr)
{
    using namespace softfloat::internals;
    uint32_t* const zWPtr = reinterpret_cast<uint32_t*>(zPtr);
    uint16_t const uiA64 = aPtr->signExp;
    bool const sign = signExtF80UI64(uiA64);
    int32_t exp = expExtF80UI64(uiA64);
    uint64_t sig = aPtr->signif;
    zWPtr[indexWord(4, 0)] = 0;

    if (exp == 0x7FFF) {
        if (sig & UINT64_C(0x7FFFFFFFFFFFFFFF)) {
            softfloat_commonNaNToF128M(softfloat_extF80MToCommonNaN(*aPtr), zWPtr);
            return;
        }

        zWPtr[indexWord(4, 3)] = packToF128UI96(sign, 0x7FFF, 0);
        zWPtr[indexWord(4, 2)] = 0;
        zWPtr[indexWord(4, 1)] = 0;
        return;
    }

    if (exp) {
        --exp;
    }

    if (0 == (sig & UINT64_C(0x8000000000000000))) {
        if (0 == sig) {
            zWPtr[indexWord(4, 3)] = packToF128UI96(sign, 0, 0);
            zWPtr[indexWord(4, 2)] = 0;
            zWPtr[indexWord(4, 1)] = 0;
            return;
        }

        exp += softfloat_normExtF80SigM(&sig);
    }

    zWPtr[indexWord(4, 1)] = static_cast<uint32_t>(sig) << 17;
    sig >>= 15;
    zWPtr[indexWord(4, 2)] = static_cast<uint32_t>(sig);

    if (exp < 0) {
        zWPtr[indexWordHi(4)] = sig >> 32;
        softfloat_shiftRight96M(
            &zWPtr[indexMultiwordHi(4, 3)],
            static_cast<uint8_t>(-exp),
            &zWPtr[indexMultiwordHi(4, 3)]);
        exp = 0;
        sig = static_cast<uint64_t>(zWPtr[indexWordHi(4)]) << 32;
    }

    zWPtr[indexWordHi(4)] = packToF128UI96(sign, static_cast<unsigned>(exp), sig >> 32);
    return;
}
