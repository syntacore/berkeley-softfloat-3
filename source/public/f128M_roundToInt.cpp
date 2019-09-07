
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

/**
@todo Make roundingMode of softfloat_round_mode type
*/
void
f128M_roundToInt(float128_t const* const aPtr,
                 uint8_t const roundingMode,
                 bool const exact,
                 float128_t* const zPtr)
{
#ifdef SOFTFLOAT_FAST_INT64
    *zPtr = f128_roundToInt(*aPtr, roundingMode, exact);
#else
    using namespace softfloat::internals::slow_int64;

    typedef uint32_t el_type;
    auto const aWPtr = reinterpret_cast<el_type const*>(aPtr);
    auto const zWPtr = reinterpret_cast<el_type*>(zPtr);
    auto const ui96 = aWPtr[indexWordHi(4)];
    auto const exp = expF128UI96(ui96);

    if (exp < 0x3FFF) {
        zWPtr[indexWord(4, 2)] = 0;
        zWPtr[indexWord(4, 1)] = 0;
        zWPtr[indexWord(4, 0)] = 0;

        auto const tmp = aWPtr[indexWord(4, 2)];
        auto const sigExtra = 0 == tmp ? aWPtr[indexWord(4, 1)] | aWPtr[indexWord(4, 0)] : tmp;

        if (0 == sigExtra && 0 == (ui96 & 0x7FFFFFFF)) {
            zWPtr[indexWordHi(4)] = ui96;
            return;
        }

        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        bool const sign = is_sign(ui96);

        switch (roundingMode) {
        case softfloat_round_near_even:
            if (0 == fracF128UI96(ui96) && 0 == sigExtra) {
                break;
            }

            [[fallthrough]];

        case softfloat_round_near_maxMag:
            if (0x3FFE == exp) {
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
        if (0x7FFF == exp && 0 != (fracF128UI96(ui96) || (aWPtr[indexWord(4, 2)] | aWPtr[indexWord(4, 1)] | aWPtr[indexWord(4, 0)]))) {
            softfloat_propagateNaNF128M(aWPtr, 0, zWPtr);
            return;
        }

        zWPtr[indexWord(4, 2)] = aWPtr[indexWord(4, 2)];
        zWPtr[indexWord(4, 1)] = aWPtr[indexWord(4, 1)];
        zWPtr[indexWord(4, 0)] = aWPtr[indexWord(4, 0)];
        zWPtr[indexWordHi(4)] = ui96;
        return;
    }

    /**
    @todo Warning   C4244   '=': conversion from 'int32_t' to 'uint8_t', possible loss of data
    */
    bool const roundNear =
        softfloat_round_near_maxMag == roundingMode ||
        softfloat_round_near_even == roundingMode;

    auto bitPos = 0x406F - exp - !!roundNear;
    auto index = indexWordLo(4);
    auto const lastIndex = indexWordHi(4);
    bool extra = false;
    el_type wordA;

    for (;;) {
        wordA = aWPtr[index];

        if (bitPos < CHAR_BIT * sizeof(el_type)) {
            break;
        }

        if (0 != wordA) {
            extra = true;
        }

        zWPtr[index] = 0;
        index += wordIncr;
        bitPos -= CHAR_BIT * sizeof(el_type);
    }

    assert(0 <= bitPos && bitPos < CHAR_BIT * sizeof(el_type));

    el_type bit = el_type(1) << bitPos;
    bool carry;
    el_type extrasMask;
    el_type wordZ;

    if (roundNear) {
        wordZ = wordA + bit;
        carry = wordZ < wordA;
        bit <<= 1;
        extrasMask = bit - 1;

        if (
            softfloat_round_near_even == roundingMode &&
            !extra &&
            0 == (wordZ & extrasMask)
        ) {
            if (0 == bit) {
                zWPtr[index] = wordZ;
                index += wordIncr;
                auto const wordZ_2 = aWPtr[index] + !!carry;
                carry = carry && 0 == wordZ_2;
                zWPtr[index] = wordZ_2 & ~1;

                while (index != lastIndex) {
                    index += wordIncr;
                    auto const wordZ_1 = aWPtr[index] + !!carry;
                    zWPtr[index] = wordZ_1;
                    carry = carry && 0 == wordZ_1;
                }

                if (exact && 0 != softfloat_compare128M(aWPtr, zWPtr)) {
                    softfloat_raiseFlags(softfloat_flag_inexact);
                }

                return;
            }

            wordZ &= ~bit;
        }
    } else {
        extrasMask = bit - 1;
        wordZ = wordA;
        carry = false;

        if (softfloat_round_minMag != roundingMode && is_sign(ui96) != (softfloat_round_max == roundingMode)) {
            if (extra || 0 != (wordA & extrasMask)) {
                wordZ += bit;
                carry = wordZ < wordA;
            }
        }
    }

    zWPtr[index] = wordZ & ~extrasMask;

    while (index != lastIndex) {
        index += wordIncr;
        auto const wordZ_1 = aWPtr[index] + !!carry;
        zWPtr[index] = wordZ_1;
        carry = carry && 0 == wordZ_1;
    }

    if (exact && 0 != softfloat_compare128M(aWPtr, zWPtr)) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }
#endif
}
