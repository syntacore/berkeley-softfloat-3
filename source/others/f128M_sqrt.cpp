
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

void f128M_sqrt(const float128_t* aPtr, float128_t* zPtr)
{

    *zPtr = f128_sqrt(*aPtr);

}

#else

void
f128M_sqrt(float128_t const* const aPtr,
           float128_t* const zPtr)
{
    uint32_t const* aWPtr = reinterpret_cast<uint32_t const*>(aPtr);
    uint32_t* const zWPtr = reinterpret_cast<uint32_t*>(zPtr);
    uint32_t const uiA96 = aWPtr[indexWordHi(4)];
    bool const signA = signF128UI96(uiA96);
    int32_t const rawExpA = expF128UI96(uiA96);

    if (0x7FFF == rawExpA) {
        if (fracF128UI96(uiA96) || 0 != (aWPtr[indexWord(4, 2)] | aWPtr[indexWord(4, 1)] | aWPtr[indexWord(4, 0)])) {
            softfloat_propagateNaNF128M(aWPtr, 0, zWPtr);
            return;
        }

        if (!signA) {
            zWPtr[indexWordHi(4)] = uiA96;
            zWPtr[indexWord(4, 2)] = aWPtr[indexWord(4, 2)];
            zWPtr[indexWord(4, 1)] = aWPtr[indexWord(4, 1)];
            zWPtr[indexWord(4, 0)] = aWPtr[indexWord(4, 0)];
            return;
        }

        softfloat_invalidF128M(zWPtr);
        return;
    }

    uint32_t rem[6];
    int32_t expA = softfloat_shiftNormSigF128M(aWPtr, static_cast<uint8_t>(13 - (rawExpA & 1)), rem);

    if (-128 == expA) {
        zWPtr[indexWordHi(4)] = uiA96;
        zWPtr[indexWord(4, 2)] = aWPtr[indexWord(4, 2)];
        zWPtr[indexWord(4, 1)] = aWPtr[indexWord(4, 1)];
        zWPtr[indexWord(4, 0)] = aWPtr[indexWord(4, 0)];
        return;
    }

    if (signA) {
        softfloat_invalidF128M(zWPtr);
        return;
    }

    /*
    `sig32Z' is guaranteed to be a lower bound on the square root of
    `sig32A', which makes `sig32Z' also a lower bound on the square root of
    `sigA'.
    */
    int32_t const expZ = ((expA - 0x3FFF) >> 1) + 0x3FFE;
    expA &= 1;
    uint64_t rem64 = static_cast<uint64_t>(rem[indexWord(4, 3)]) << 32 | rem[indexWord(4, 2)];
    uint32_t sig32A;

    if (expA) {
        if (!rawExpA) {
            softfloat_shortShiftRight128M(rem, 1, rem);
            rem64 >>= 1;
        }

        sig32A = static_cast<uint32_t>(rem64 >> 29);
    } else {
        sig32A = static_cast<uint32_t>(rem64 >> 30);
    }

    uint32_t const recipSqrt32 = softfloat_approxRecipSqrt32_1(expA, sig32A);
    uint32_t sig32Z = (static_cast<uint64_t>(sig32A) * recipSqrt32) >> 32;

    if (expA) {
        sig32Z >>= 1;
    }

    uint32_t qs[3];
    qs[2] = sig32Z;
    rem64 -= (uint64_t)sig32Z * sig32Z;
    rem[indexWord(4, 3)] = rem64 >> 32;
    rem[indexWord(4, 2)] = (uint32_t)rem64;

    uint32_t q = (static_cast<uint32_t>(rem64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    qs[1] = q;
    uint64_t sig64Z = (static_cast<uint64_t>(sig32Z) << 32) + (static_cast<uint64_t>(q) << 3);
    uint64_t x64 = (static_cast<uint64_t>(sig32Z) << 32) + sig64Z;
    uint32_t term[5];
    term[indexWord(4, 3)] = 0;
    term[indexWord(4, 2)] = x64 >> 32;
    term[indexWord(4, 1)] = (uint32_t)x64;
    term[indexWord(4, 0)] = 0;
    uint32_t y[5];
    softfloat_remStep128MBy32(rem, 29, term, q, y);
    rem64 = (uint64_t)y[indexWord(4, 3)] << 32 | y[indexWord(4, 2)];

    q = (static_cast<uint32_t>(rem64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    sig64Z <<= 1;

    /* Repeating this loop is a rare occurrence. */
    uint32_t rem32;

    for (;;) {
        x64 = sig64Z + (q >> 26);
        term[indexWord(4, 2)] = x64 >> 32;
        term[indexWord(4, 1)] = (uint32_t)x64;
        term[indexWord(4, 0)] = q << 6;
        term[indexWord(4, 3)] = 0;
        softfloat_remStep128MBy32(y, 29, term, q, &rem[indexMultiwordHi(6, 4)]);
        rem32 = rem[indexWordHi(6)];

        if (!(rem32 & 0x80000000)) {
            break;
        }

        --q;
    }

    qs[0] = q;
    rem64 = (uint64_t)rem32 << 32 | rem[indexWord(6, 4)];

    q = ((static_cast<uint32_t>(rem64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32) + 2;
    x64 = static_cast<uint64_t>(q) << 27;
    y[indexWord(5, 0)] = static_cast<uint32_t>(x64);
    x64 = ((uint64_t)qs[0] << 24) + (x64 >> 32);
    y[indexWord(5, 1)] = static_cast<uint32_t>(x64);
    x64 = (static_cast<uint64_t>(qs[1]) << 21) + (x64 >> 32);
    y[indexWord(5, 2)] = static_cast<uint32_t>(x64);
    x64 = ((uint64_t)qs[2] << 18) + (x64 >> 32);
    y[indexWord(5, 3)] = static_cast<uint32_t>(x64);
    y[indexWord(5, 4)] = x64 >> 32;

    if ((q & 0xF) <= 2) {
        q &= ~3;
        y[indexWordLo(5)] = q << 27;
        term[indexWord(5, 4)] = 0;
        term[indexWord(5, 3)] = 0;
        term[indexWord(5, 2)] = 0;
        term[indexWord(5, 1)] = q >> 6;
        term[indexWord(5, 0)] = q << 26;
        softfloat_sub160M(y, term, term);
        rem[indexWord(6, 1)] = 0;
        rem[indexWord(6, 0)] = 0;
        softfloat_remStep160MBy32(
            &rem[indexMultiwordLo(6, 5)],
            14,
            term,
            q,
            &rem[indexMultiwordLo(6, 5)]
        );
        rem32 = rem[indexWord(6, 4)];

        if (rem32 & 0x80000000) {
            softfloat_sub1X160M(y);
        } else if (0 != rem32 || 0 != rem[indexWord(6, 0)] || 0 != rem[indexWord(6, 1)] || 0 != (rem[indexWord(6, 3)] | rem[indexWord(6, 2)])) {
            y[indexWordLo(5)] |= 1;
        }
    }

    softfloat_roundPackMToF128M(0, expZ, y, zWPtr);
    return;
}

#endif

