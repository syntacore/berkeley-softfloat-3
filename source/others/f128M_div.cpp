
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
f128M_div(const float128_t* aPtr, const float128_t* bPtr, float128_t* zPtr)
{

    *zPtr = f128_div(*aPtr, *bPtr);

}

#else

void
f128M_div(float128_t const* const aPtr,
          float128_t const* const bPtr,
          float128_t* const zPtr)
{
    using namespace softfloat;
    uint32_t const* const aWPtr = reinterpret_cast<uint32_t const*>(aPtr);
    uint32_t const* const bWPtr = reinterpret_cast<uint32_t const*>(bPtr);
    uint32_t* const zWPtr = reinterpret_cast<uint32_t*>(zPtr);

    uint32_t const uiA96 = aWPtr[indexWordHi(4)];
    bool const signA = signF128UI96(uiA96);
    int32_t expA = expF128UI96(uiA96);
    uint32_t const uiB96 = bWPtr[indexWordHi(4)];
    bool const signB = signF128UI96(uiB96);
    int32_t expB = expF128UI96(uiB96);
    bool const signZ = signA != signB;

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (softfloat_tryPropagateNaNF128M(aWPtr, bWPtr, zWPtr)) {
            return;
        }

        if (expA == 0x7FFF) {
            if (expB == 0x7FFF) {
                softfloat_invalidF128M(zWPtr);
                return;
            }

            zWPtr[indexWordHi(4)] = packToF128UI96(signZ, 0x7FFF, 0);
            zWPtr[indexWord(4, 2)] = 0;
            zWPtr[indexWord(4, 1)] = 0;
            zWPtr[indexWord(4, 0)] = 0;
            return;
        }

        zWPtr[indexWordHi(4)] = packToF128UI96(signZ, 0, 0);
        zWPtr[indexWord(4, 2)] = 0;
        zWPtr[indexWord(4, 1)] = 0;
        zWPtr[indexWord(4, 0)] = 0;
        return;
    }

    uint32_t y[5];
    expA = softfloat_shiftNormSigF128M(aWPtr, 13, y);
    uint32_t sigB[4];
    expB = softfloat_shiftNormSigF128M(bWPtr, 13, sigB);

    if (expA == -128) {
        if (expB == -128) {
            softfloat_invalidF128M(zWPtr);
            return;
        }

        zWPtr[indexWordHi(4)] = packToF128UI96(signZ, 0, 0);
        zWPtr[indexWord(4, 2)] = 0;
        zWPtr[indexWord(4, 1)] = 0;
        zWPtr[indexWord(4, 0)] = 0;
        return;
    }

    if (expB == -128) {
        softfloat_raiseFlags(softfloat_flag_infinite);
        zWPtr[indexWordHi(4)] = packToF128UI96(signZ, 0x7FFF, 0);
        zWPtr[indexWord(4, 2)] = 0;
        zWPtr[indexWord(4, 1)] = 0;
        zWPtr[indexWord(4, 0)] = 0;
        return;
    }

    int32_t expZ = expA - expB + 0x3FFE;

    if (softfloat_compare128M(y, sigB) < 0) {
        --expZ;
        softfloat_add128M(y, y, y);
    }

    uint32_t const recip32 = softfloat_approxRecip32_1((static_cast<uint64_t>(sigB[indexWord(4, 3)]) << 32 | sigB[indexWord(4, 2)]) >> 30);
    int ix = 3;
    uint64_t q64;
    uint32_t qs[3];
    uint32_t q;

    for (;;) {
        q64 = (uint64_t)y[indexWordHi(4)] * recip32;
        q = (q64 + 0x80000000) >> 32;
        --ix;

        if (ix < 0) {
            break;
        }

        softfloat_remStep128MBy32(y, 29, sigB, q, y);

        if (y[indexWordHi(4)] & 0x80000000) {
            --q;
            softfloat_add128M(y, sigB, y);
        }

        qs[ix] = q;
    }

    if (((q + 1) & 7) < 2) {
        softfloat_remStep128MBy32(y, 29, sigB, q, y);

        if (y[indexWordHi(4)] & 0x80000000) {
            --q;
            softfloat_add128M(y, sigB, y);
        } else if (softfloat_compare128M(sigB, y) <= 0) {
            ++q;
            softfloat_sub128M(y, sigB, y);
        }

        if (y[indexWordLo(4)] || y[indexWord(4, 1)] || (y[indexWord(4, 2)] | y[indexWord(4, 3)])) {
            q |= 1;
        }
    }

    q64 = (uint64_t)q << 28;
    y[indexWord(5, 0)] = (uint32_t)q64;
    q64 = ((uint64_t)qs[0] << 25) + (q64 >> 32);
    y[indexWord(5, 1)] = (uint32_t)q64;
    q64 = ((uint64_t)qs[1] << 22) + (q64 >> 32);
    y[indexWord(5, 2)] = (uint32_t)q64;
    q64 = ((uint64_t)qs[2] << 19) + (q64 >> 32);
    y[indexWord(5, 3)] = (uint32_t)q64;
    y[indexWord(5, 4)] = q64 >> 32;
    softfloat_roundPackMToF128M(signZ, expZ, y, zWPtr);
    return;
}

#endif
