
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015 The Regents of the University of
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

#include "internals.h"

#include "softfloat/functions.h"

void
softfloat_roundPackMToExtF80M(
    bool sign,
    int32_t exp,
    uint32_t *extSigPtr,
    uint8_t roundingPrecision,
   /** @bug use extFloat80_t */
   struct extFloat80M *zSPtr
)
{
    uint64_t roundIncrement;
    uint64_t roundMask;
    uint64_t roundBits;
    uint32_t sigExtra;
    bool doIncrement;


    uint8_t const roundingMode = softfloat_roundingMode;
    bool const roundNearEven = roundingMode == softfloat_round_near_even;
    uint64_t sig =
        (uint64_t)extSigPtr[indexWord(3, 2)] << 32 |
        extSigPtr[indexWord(3, 1)];
    if (roundingPrecision == 80) {
        goto precision80;
    }
    if (roundingPrecision == 64) {
        roundIncrement = UINT64_C(0x0000000000000400);
        roundMask = UINT64_C(0x00000000000007FF);
    } else if (roundingPrecision == 32) {
        roundIncrement = UINT64_C(0x0000008000000000);
        roundMask = UINT64_C(0x000000FFFFFFFFFF);
    } else {
        goto precision80;
    }

    if (extSigPtr[indexWordLo(3)]) {
        sig |= 1;
    }
    if (!roundNearEven && (roundingMode != softfloat_round_near_maxMag)) {
        roundIncrement =
            roundingMode == (sign ? softfloat_round_min : softfloat_round_max) ?
            roundMask : 0;
    }
    roundBits = sig & roundMask;

    if (0x7FFD <= (uint32_t)(exp - 1)) {
        if (exp <= 0) {
            bool const isTiny =
                softfloat_detectTininess == softfloat_tininess_beforeRounding ||
                exp < 0 ||
                sig <= (uint64_t)(sig + roundIncrement);
            sig = softfloat_shiftRightJam64(sig, 1 - exp);
            roundBits = sig & roundMask;
            if (roundBits) {
                if (isTiny) {
                    softfloat_raiseFlags(softfloat_flag_underflow);
                }
                softfloat_raiseFlags(softfloat_flag_inexact);
            }
            sig += roundIncrement;
            exp = 0 != (sig & INT64_MIN);
            roundIncrement = roundMask + 1;
            if (roundNearEven && (roundBits << 1 == roundIncrement)) {
                roundMask |= roundIncrement;
            }
            sig &= ~roundMask;
            goto packReturn;
        }
        if (0x7FFE < exp || (exp == 0x7FFE && (uint64_t)(sig + roundIncrement) < sig)) {
            goto overflow;
        }
    }

    if (roundBits) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }
    sig += roundIncrement;
    if (sig < roundIncrement) {
        ++exp;
        sig = (uint64_t)INT64_MIN;
    }
    roundIncrement = roundMask + 1;
    if (roundNearEven && (roundBits << 1 == roundIncrement)) {
        roundMask |= roundIncrement;
    }
    sig &= ~roundMask;
    goto packReturn;

precision80:
    sigExtra = extSigPtr[indexWordLo(3)];
    doIncrement = (0x80000000 <= sigExtra);
    if (!roundNearEven && (roundingMode != softfloat_round_near_maxMag)) {
        doIncrement =
            roundingMode == (sign ? softfloat_round_min : softfloat_round_max) &&
            sigExtra;
    }

    if (0x7FFD <= (uint32_t)(exp - 1)) {
        if (exp <= 0) {
            bool const isTiny =
                softfloat_detectTininess == softfloat_tininess_beforeRounding ||
                exp < 0 ||
                !doIncrement ||
                sig < UINT64_MAX;
            softfloat_shiftRightJam96M(extSigPtr, 1 - exp, extSigPtr);
            sigExtra = extSigPtr[indexWordLo(3)];
            if (sigExtra) {
                if (isTiny) {
                    softfloat_raiseFlags(softfloat_flag_underflow);
                }
                softfloat_raiseFlags(softfloat_flag_inexact);
            }
            doIncrement = (0x80000000 <= sigExtra);
            if (!roundNearEven && roundingMode != softfloat_round_near_maxMag) {
                doIncrement =
                    roundingMode == (sign ? softfloat_round_min : softfloat_round_max) &&
                    sigExtra;
            }
            exp = 0;
            sig =
                (uint64_t)extSigPtr[indexWord(3, 2)] << 32 |
                extSigPtr[indexWord(3, 1)];
            if (doIncrement) {
                ++sig;
                sig &= ~(uint64_t)(!(sigExtra & 0x7FFFFFFF) & roundNearEven);
                exp = ((sig & UINT64_C(0x8000000000000000)) != 0);
            }
            goto packReturn;
        }

        if (0x7FFE < exp || (exp == 0x7FFE && sig == UINT64_MAX && doIncrement)) {
            roundMask = 0;
overflow:
            softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);
            if (roundNearEven ||
                roundingMode == softfloat_round_near_maxMag ||
                roundingMode == (sign ? softfloat_round_min : softfloat_round_max)
            ) {
                exp = INT16_MAX;
                sig = (uint64_t)INT64_MIN;
            } else {
                exp = INT16_MAX - 1;
                sig = ~roundMask;
            }
            goto packReturn;
        }
    }

    if (sigExtra) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }
    if (doIncrement) {
        ++sig;
        if (!sig) {
            ++exp;
            sig = (uint64_t)INT64_MIN;
        } else {
            sig &= ~(uint64_t)(!(sigExtra & INT32_MAX) & roundNearEven);
        }
    }

packReturn:
    /** @todo Warning	C4244	'=': conversion from 'int32_t' to 'uint16_t', possible loss of data */
    zSPtr->signExp = packToExtF80UI64(sign, exp);
    zSPtr->signif = sig;
}
