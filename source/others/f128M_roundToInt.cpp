
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

#include "softfloat/functions.h"

#include "internals.hpp"
#include "specialize.hpp"

/** @todo split to different implementations */
#ifdef SOFTFLOAT_FAST_INT64

void
f128M_roundToInt(float128_t const* aPtr,
                 uint8_t roundingMode,
                 bool exact,
                 float128_t* zPtr)
{
    *zPtr = f128_roundToInt(*aPtr, roundingMode, exact);
}

#else

void
f128M_roundToInt(float128_t const* aPtr,
                 uint8_t roundingMode,
                 bool exact,
                 float128_t* zPtr)
{
    using namespace softfloat;
    uint32_t sigExtra;
    bool sign;
    unsigned int index, lastIndex;
    bool extra;
    uint32_t wordA, wordZ;
    uint32_t extrasMask;


    uint32_t const* const aWPtr = reinterpret_cast<uint32_t const*>(aPtr);
    uint32_t* const zWPtr = reinterpret_cast<uint32_t*>(zPtr);

    uint32_t const ui96 = aWPtr[indexWordHi(4)];
    auto const exp = expF128UI96(ui96);

    if (exp < 0x3FFF) {
        zWPtr[indexWord(4, 2)] = 0;
        zWPtr[indexWord(4, 1)] = 0;
        zWPtr[indexWord(4, 0)] = 0;
        sigExtra = aWPtr[indexWord(4, 2)];

        if (!sigExtra) {
            sigExtra = aWPtr[indexWord(4, 1)] | aWPtr[indexWord(4, 0)];
        }

        if (!sigExtra && !(ui96 & 0x7FFFFFFF)) {
            zWPtr[indexWordHi(4)] = ui96;
            return;
        }

        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        sign = signF128UI96(ui96);

        switch (roundingMode) {
        case softfloat_round_near_even:
            if (!fracF128UI96(ui96) && !sigExtra) {
                break;
            }

        case softfloat_round_near_maxMag:
            if (exp == 0x3FFE) {
                zWPtr[indexWordHi(4)] = packToF128UI96(sign, 0x3FFF, 0);
                return;
            }

            break;

        case softfloat_round_min:
            if (sign) {
                zWPtr[indexWordHi(4)] = packToF128UI96(sign, 0x3FFF, 0);
                return;
            }

            break;

        case softfloat_round_max:
            if (!sign) {
                zWPtr[indexWordHi(4)] = packToF128UI96(sign, 0x3FFF, 0);
                return;
            }

            break;
        }

        zWPtr[indexWordHi(4)] = packToF128UI96(sign, 0, 0);
        return;
    }

    if (0x406F <= exp) {
        if (0x7FFF == exp && (fracF128UI96(ui96) || (aWPtr[indexWord(4, 2)] | aWPtr[indexWord(4, 1)] | aWPtr[indexWord(4, 0)]))) {
            softfloat_propagateNaNF128M(aWPtr, 0, zWPtr);
            return;
        }

        zWPtr[indexWord(4, 2)] = aWPtr[indexWord(4, 2)];
        zWPtr[indexWord(4, 1)] = aWPtr[indexWord(4, 1)];
        zWPtr[indexWord(4, 0)] = aWPtr[indexWord(4, 0)];
        zWPtr[indexWordHi(4)] = ui96;
        return;
    }

    /** @todo Warning   C4244   '=': conversion from 'int32_t' to 'uint8_t', possible loss of data */
    bool const roundNear =
        softfloat_round_near_maxMag == roundingMode ||
        softfloat_round_near_even == roundingMode;
    auto bitPos = 0x406F - exp - !!roundNear;
    index = indexWordLo(4);
    lastIndex = indexWordHi(4);
    extra = 0;

    for (;;) {
        wordA = aWPtr[index];

        if (bitPos < 32) {
            break;
        }

        if (wordA) {
            extra = 1;
        }

        zWPtr[index] = 0;
        index += wordIncr;
        bitPos -= 32;
    }

    uint32_t bit = UINT32_C(1) << bitPos;

    bool carry;

    if (roundNear) {
        wordZ = wordA + bit;
        carry = wordZ < wordA;
        bit <<= 1;
        extrasMask = bit - 1;

        if (softfloat_round_near_even == roundingMode &&
                !extra &&
                0 == (wordZ & extrasMask)) {
            if (!bit) {
                zWPtr[index] = wordZ;
                index += wordIncr;
                wordZ = aWPtr[index] + !!carry;
                carry = carry && !wordZ;
                zWPtr[index] = wordZ & ~1;
                goto propagateCarry;
            }

            wordZ &= ~bit;
        }
    } else {
        extrasMask = bit - 1;
        wordZ = wordA;
        carry = false;

        if (softfloat_round_minMag != roundingMode && (signF128UI96(ui96) ^ (roundingMode == softfloat_round_max))) {
            if (extra || (wordA & extrasMask)) {
                wordZ += bit;
                carry = wordZ < wordA;
            }
        }
    }

    wordZ &= ~extrasMask;
    zWPtr[index] = wordZ;

propagateCarry:

    while (index != lastIndex) {
        index += wordIncr;
        wordZ = aWPtr[index] + !!carry;
        zWPtr[index] = wordZ;
        carry = carry && !wordZ;
    }

    if (exact && (softfloat_compare128M(aWPtr, zWPtr) != 0)) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }
}

#endif
