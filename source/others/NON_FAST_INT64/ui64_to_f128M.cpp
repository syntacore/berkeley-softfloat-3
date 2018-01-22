
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All Rights Reserved.

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

void
ui64_to_f128M(uint64_t a,
              float128_t* zPtr)
{
    using namespace softfloat::internals;
    uint32_t* const zWPtr = (uint32_t*)zPtr;
    uint32_t uiZ96 = 0;
    uint32_t uiZ64 = 0;
    zWPtr[indexWord(4, 1)] = 0;
    zWPtr[indexWord(4, 0)] = 0;

    if (a) {
        uint8_t const shiftDist = softfloat_countLeadingZeros64(a) + 17u;

        if (shiftDist < 32) {
            uint32_t* const ptr = zWPtr + indexMultiwordHi(4, 3);
            ptr[indexWord(3, 2)] = 0;
            ptr[indexWord(3, 1)] = a >> 32;
            ptr[indexWord(3, 0)] = static_cast<uint32_t>(a);
            softfloat_shortShiftLeft96M(ptr, shiftDist, ptr);
            ptr[indexWordHi(3)] = packToF128UI96(0, 0x404Eu - shiftDist, ptr[indexWordHi(3)]);
            return;
        }

        a <<= shiftDist - 32;
        uiZ96 = packToF128UI96(0, 0x404Eu - shiftDist, a >> 32);
        uiZ64 = static_cast<uint32_t>(a);
    }

    zWPtr[indexWord(4, 3)] = uiZ96;
    zWPtr[indexWord(4, 2)] = uiZ64;
}
