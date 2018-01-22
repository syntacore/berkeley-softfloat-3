
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

void
f128M_add(const float128_t* aPtr,
          const float128_t* bPtr,
          float128_t* zPtr)
{
    using namespace softfloat::internals;
    const uint64_t* aWPtr, *bWPtr;
    uint64_t uiA64, uiA0;
    bool signA;
    uint64_t uiB64, uiB0;
    bool signB;
    aWPtr = (const uint64_t*)aPtr;
    bWPtr = (const uint64_t*)bPtr;
    uiA64 = aWPtr[indexWord(2, 1)];
    uiA0 = aWPtr[indexWord(2, 0)];
    signA = signF128UI64(uiA64);
    uiB64 = bWPtr[indexWord(2, 1)];
    uiB0 = bWPtr[indexWord(2, 0)];
    signB = signF128UI64(uiB64);

    if (signA == signB) {
        *zPtr = softfloat_addMagsF128(uiA64, uiA0, uiB64, uiB0, signA);
    } else {
        *zPtr = softfloat_subMagsF128(uiA64, uiA0, uiB64, uiB0, signA);
    }
}
