
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

bool
f128M_eq(float128_t const* aPtr,
         float128_t const* bPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    return f128_eq(*aPtr, *bPtr);
#else
    using namespace softfloat::internals::slow_int64;

    uint32_t const* const aWPtr = reinterpret_cast<uint32_t const*>(aPtr);
    uint32_t const* const bWPtr = reinterpret_cast<uint32_t const*>(bPtr);
    uint32_t wordA = aWPtr[indexWord(4, 2)];
    uint32_t wordB = bWPtr[indexWord(4, 2)];

    if (wordA == wordB) {
        uint32_t const uiA96 = aWPtr[indexWordHi(4)];
        uint32_t const uiB96 = bWPtr[indexWordHi(4)];
        bool possibleOppositeZeros = false;

        if (uiA96 != uiB96) {
            possibleOppositeZeros = 0 == ((uiA96 | uiB96) & 0x7FFFFFFF);

            if (!possibleOppositeZeros) {
                if (f128M_isSignalingNaN(aPtr) || f128M_isSignalingNaN(bPtr)) {
                    softfloat_raiseFlags(softfloat_flag_invalid);
                }

                return false;
            }
        }

        auto const mashWord = wordA | wordB;
        auto const wordA_2 = aWPtr[indexWord(4, 1)];
        auto const wordB_2 = bWPtr[indexWord(4, 1)];

        if (wordA_2 == wordB_2) {
            auto const mashWord_1 = mashWord | wordA_2 | wordB_2;
            auto const wordA_1 = aWPtr[indexWord(4, 0)];
            auto const wordB_1 = bWPtr[indexWord(4, 0)];

            if (
                wordA_1 == wordB_1 &&
                !(possibleOppositeZeros && 0 != (mashWord_1 | wordA_1 | wordB_1)) &&
                !softfloat_isNaNF128M(aWPtr) && !softfloat_isNaNF128M(bWPtr)
            ) {
                return true;
            }
        }
    }

    if (f128M_isSignalingNaN(aPtr) || f128M_isSignalingNaN(bPtr)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return false;
#endif
}
