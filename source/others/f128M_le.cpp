
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
f128M_le(float128_t const *const aPtr,
         float128_t const *const bPtr)
{
#ifdef SOFTFLOAT_FAST_INT64
    return f128_le(*aPtr, *bPtr);
#else
    using namespace softfloat::internals::slow_int64;

    uint32_t const* aWPtr = reinterpret_cast<uint32_t const *>(aPtr);
    uint32_t const* bWPtr = reinterpret_cast<uint32_t const *>(bPtr);

    if (softfloat_isNaNF128M(aWPtr) || softfloat_isNaNF128M(bWPtr)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return false;
    }

    uint32_t const uiA96 = aWPtr[indexWordHi(4)];
    uint32_t const uiB96 = bWPtr[indexWordHi(4)];
    bool const signA = is_sign(uiA96);
    bool const signB = is_sign(uiB96);

    if (signA != signB) {
        if (signA) {
            return true;
        }

        if (0 != ((uiA96 | uiB96) & 0x7FFFFFFF)) {
            return false;
        }

        uint32_t wordA = aWPtr[indexWord(4, 2)];
        uint32_t wordB = bWPtr[indexWord(4, 2)];

        if (0 != (wordA | wordB)) {
            return false;
        }

        wordA = aWPtr[indexWord(4, 1)];
        wordB = bWPtr[indexWord(4, 1)];

        if (0 != (wordA | wordB)) {
            return false;
        }

        wordA = aWPtr[indexWord(4, 0)];
        wordB = bWPtr[indexWord(4, 0)];
        return 0 == (wordA | wordB);
    }

    if (signA) {
        std::swap(aWPtr, bWPtr);
    }

    return (softfloat_compare128M(aWPtr, bWPtr) <= 0);
#endif
}
