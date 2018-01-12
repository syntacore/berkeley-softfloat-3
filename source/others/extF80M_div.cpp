
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All rights reserved.
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
#include "specialize.hpp"

/** @todo split to different implementations */
#ifdef SOFTFLOAT_FAST_INT64

void
extF80M_div(const extFloat80_t *aPtr, const extFloat80_t *bPtr, extFloat80_t *zPtr)
{
    *zPtr = extF80_div(*aPtr, *bPtr);
}

#else

void
extF80M_div(extFloat80_t const *aSPtr, extFloat80_t const *bSPtr, extFloat80_t *zSPtr)
{
    uint16_t const uiA64 = aSPtr->signExp;
    int32_t expA = expExtF80UI64(uiA64);
    uint16_t const uiB64 = bSPtr->signExp;
    int32_t expB = expExtF80UI64(uiB64);
    bool const signZ = signExtF80UI64(uiA64) ^ signExtF80UI64(uiB64);
    if (expA == INT16_MAX || expB == INT16_MAX) {
        if (softfloat_tryPropagateNaNExtF80M(aSPtr, bSPtr, zSPtr)) {
            return;
        } else if (expA != INT16_MAX) {
            zSPtr->signExp = packToExtF80UI64(signZ, 0);
            zSPtr->signif = 0;
            return;
        } else if (expB == INT16_MAX) {
            softfloat_invalidExtF80M(zSPtr);
            return;
        } else {
            zSPtr->signExp = packToExtF80UI64(signZ, 0x7FFF);
            zSPtr->signif = UINT64_C(0x8000000000000000);
            return;
        }
    } else {
        uint64_t sigA = aSPtr->signif;
        uint64_t x64 = bSPtr->signif;
        if (!expB) {
            expB = 1;
        }
        if (!(x64 & UINT64_C(0x8000000000000000))) {
            if (!x64) {
                if (!sigA) {
                    softfloat_invalidExtF80M(zSPtr);
                    return;
                } else {
                    softfloat_raiseFlags(softfloat_flag_infinite);
                    zSPtr->signExp = packToExtF80UI64(signZ, 0x7FFF);
                    zSPtr->signif = UINT64_C(0x8000000000000000);
                    return;
                }
            } else {
                expB += softfloat_normExtF80SigM(&x64);
            }
        }
        if (!expA) {
            expA = 1;
        }
        if (!(sigA & UINT64_C(0x8000000000000000))) {
            if (!sigA) {
                zSPtr->signExp = packToExtF80UI64(signZ, 0);
                zSPtr->signif = 0;
                return;
            } else {
                expA += softfloat_normExtF80SigM(&sigA);
            }
        }

        int32_t expZ = expA - expB + 0x3FFF;
        uint8_t shiftDist = 29;
        if (sigA < x64) {
            --expZ;
            shiftDist = 30;
        }
        uint32_t y[3];
        softfloat_shortShiftLeft64To96M(sigA, shiftDist, y);
        uint32_t const recip32 = softfloat_approxRecip32_1(x64 >> 32);
        uint32_t sigB[3];
        sigB[indexWord(3, 0)] = (uint32_t)x64 << 30;
        x64 >>= 2;
        sigB[indexWord(3, 2)] = x64 >> 32;
        sigB[indexWord(3, 1)] = (uint32_t)x64;
        int ix = 2;
        uint32_t q;
        uint32_t qs[2];
        for (;;) {
            x64 = (uint64_t)y[indexWordHi(3)] * recip32;
            q = (x64 + 0x80000000) >> 32;
            --ix;
            if (ix < 0) {
                break;
            }
            softfloat_remStep96MBy32(y, 29, sigB, q, y);
            if (y[indexWordHi(3)] & 0x80000000) {
                --q;
                softfloat_add96M(y, sigB, y);
            }
            qs[ix] = q;
        }

        if (((q + 1) & 0x3FFFFF) < 2) {
            softfloat_remStep96MBy32(y, 29, sigB, q, y);
            if (y[indexWordHi(3)] & 0x80000000) {
                --q;
                softfloat_add96M(y, sigB, y);
            } else if (softfloat_compare96M(sigB, y) <= 0) {
                ++q;
                softfloat_sub96M(y, sigB, y);
            }
            if (y[indexWordLo(3)] || y[indexWord(3, 1)] || y[indexWord(3, 2)]) {
                q |= 1;
            }
        }
        x64 = (uint64_t)q << 9;
        y[indexWord(3, 0)] = (uint32_t)x64;
        x64 = ((uint64_t)qs[0] << 6) + (x64 >> 32);
        y[indexWord(3, 1)] = (uint32_t)x64;
        y[indexWord(3, 2)] = (qs[1] << 3) + (x64 >> 32);
        softfloat_roundPackMToExtF80M(signZ, expZ, y, extF80_roundingPrecision, zSPtr);
        return;
    }
}

#endif

