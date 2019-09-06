
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
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
f128M_mul(float128_t const* const aPtr,
          float128_t const* const bPtr,
          float128_t* const zPtr)
{
#ifdef SOFTFLOAT_FAST_INT64
    *zPtr = f128_mul(*aPtr, *bPtr);
#else
    using namespace softfloat::internals;

    uint32_t const* const aWPtr = reinterpret_cast<const uint32_t*>(aPtr);
    uint32_t const* const bWPtr = reinterpret_cast<const uint32_t*>(bPtr);
    uint32_t* const zWPtr = reinterpret_cast<uint32_t*>(zPtr);

    uint32_t const uiA96 = aWPtr[indexWordHi(4)];
    int32_t expA = expF128UI96(uiA96);
    uint32_t const uiB96 = bWPtr[indexWordHi(4)];
    int32_t expB = expF128UI96(uiB96);
    bool const signZ = is_sign(uiA96) != is_sign(uiB96);

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (softfloat_tryPropagateNaNF128M(aWPtr, bWPtr, zWPtr)) {
            return;
        }

        if (0 == expA || 0 == expB) {
            uint32_t const* ptr = 0 == expA ? aWPtr : bWPtr;

            if (
                0 == fracF128UI96(ptr[indexWordHi(4)]) &&
                0 == (ptr[indexWord(4, 2)] | ptr[indexWord(4, 1)] | ptr[indexWord(4, 0)])
            ) {
                softfloat_invalidF128M(zWPtr);
                return;
            }
        }

        zWPtr[indexWordHi(4)] = packToF128UI96(signZ, 0x7FFF, 0);
        zWPtr[indexWord(4, 2)] = 0;
        zWPtr[indexWord(4, 1)] = 0;
        zWPtr[indexWord(4, 0)] = 0;
        return;
    }

    uint32_t sigA[4];

    if (0 != expA) {
        sigA[indexWordHi(4)] = fracF128UI96(uiA96) | 0x00010000;
        sigA[indexWord(4, 2)] = aWPtr[indexWord(4, 2)];
        sigA[indexWord(4, 1)] = aWPtr[indexWord(4, 1)];
        sigA[indexWord(4, 0)] = aWPtr[indexWord(4, 0)];
    } else {
        expA = softfloat_shiftNormSigF128M(aWPtr, 0, sigA);

        if (expA == -128) {
            zWPtr[indexWordHi(4)] = packToF128UI96(signZ, 0, 0);
            zWPtr[indexWord(4, 2)] = 0;
            zWPtr[indexWord(4, 1)] = 0;
            zWPtr[indexWord(4, 0)] = 0;
            return;
        }
    }

    uint32_t sigB[4];

    if (0 != expB) {
        sigB[indexWordHi(4)] = fracF128UI96(uiB96) | 0x00010000;
        sigB[indexWord(4, 2)] = bWPtr[indexWord(4, 2)];
        sigB[indexWord(4, 1)] = bWPtr[indexWord(4, 1)];
        sigB[indexWord(4, 0)] = bWPtr[indexWord(4, 0)];
    } else {
        expB = softfloat_shiftNormSigF128M(bWPtr, 0, sigB);

        if (expB == -128) {
            zWPtr[indexWordHi(4)] = packToF128UI96(signZ, 0, 0);
            zWPtr[indexWord(4, 2)] = 0;
            zWPtr[indexWord(4, 1)] = 0;
            zWPtr[indexWord(4, 0)] = 0;
            return;
        }
    }

    int32_t expZ = expA + expB - 0x4000;
    uint32_t sigProd[8];
    softfloat_mul128MTo256M(sigA, sigB, sigProd);

    if (0 != sigProd[indexWord(8, 2)] || 0 != (sigProd[indexWord(8, 1)] | sigProd[indexWord(8, 0)])) {
        sigProd[indexWord(8, 3)] |= 1;
    }

    uint32_t* const extSigZPtr = &sigProd[indexMultiwordHi(8, 5)];
    uint8_t shiftDist = 16;

    if (0 != (extSigZPtr[indexWordHi(5)] & 2)) {
        ++expZ;
        shiftDist = 15;
    }

    softfloat_shortShiftLeft160M(extSigZPtr, shiftDist, extSigZPtr);
    softfloat_roundPackMToF128M(signZ, expZ, extSigZPtr, zWPtr);
#endif
}
