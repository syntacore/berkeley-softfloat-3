
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

#include "internals.hpp"

#include "softfloat/functions.h"

namespace softfloat {

/**
@param[in] sign
@param[in] exp
@param[inout] extSigPtr
@param[in] roundingPrecision
@param[out] zSPtr
*/
void
softfloat_roundPackMToExtF80M(bool sign,
                              int32_t exp,
                              uint32_t* extSigPtr,
                              uint8_t roundingPrecision,
                              extFloat80M* zSPtr)
{
    uint64_t sig =
        static_cast<uint64_t>(extSigPtr[indexWord(3, 2)]) << 32 |
        extSigPtr[indexWord(3, 1)];

    if (64 == roundingPrecision || 32 == roundingPrecision) {
        uint64_t roundIncrement;
        uint64_t roundMask;

        if (roundingPrecision == 64) {
            roundIncrement = UINT64_C(0x0000000000000400);
            roundMask = UINT64_C(0x00000000000007FF);
        } else {
            /* roundingPrecision == 32 */
            roundIncrement = UINT64_C(0x0000008000000000);
            roundMask = UINT64_C(0x000000FFFFFFFFFF);
        }

        if (extSigPtr[indexWordLo(3)]) {
            sig |= 1;
        }

        if (softfloat_round_near_even != softfloat_roundingMode && softfloat_round_near_maxMag != softfloat_roundingMode) {
            roundIncrement =
                (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode ?
                roundMask : 0;
        }

        uint64_t roundBits = sig & roundMask;

        if (0x7FFD <= static_cast<uint32_t>(exp - 1)) {
            if (exp <= 0) {
                bool const isTiny =
                    softfloat_detectTininess == softfloat_tininess_beforeRounding ||
                    exp < 0 ||
                    sig <= static_cast<uint64_t>(sig + roundIncrement);
                sig = softfloat_shiftRightJam64(sig, static_cast<uint32_t>(1 - exp));
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

                if (softfloat_round_near_even == softfloat_roundingMode && (roundBits << 1 == roundIncrement)) {
                    roundMask |= roundIncrement;
                }

                sig &= ~roundMask;
                zSPtr->signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
                zSPtr->signif = sig;
                return;
            }

            if (0x7FFE < exp || (exp == 0x7FFE && static_cast<uint64_t>(sig + roundIncrement) < sig)) {
                softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);

                if (
                    softfloat_round_near_even == softfloat_roundingMode ||
                    softfloat_round_near_maxMag == softfloat_roundingMode ||
                    (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode
                ) {
                    exp = INT16_MAX;
                    sig = static_cast<uint64_t>(INT64_MIN);
                } else {
                    exp = INT16_MAX - 1;
                    sig = ~roundMask;
                }

                zSPtr->signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
                zSPtr->signif = sig;
                return;
            }
        }

        if (roundBits) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        sig += roundIncrement;

        if (sig < roundIncrement) {
            ++exp;
            sig = static_cast<uint64_t>(INT64_MIN);
        }

        roundIncrement = roundMask + 1;

        if (softfloat_round_near_even == softfloat_roundingMode && (roundBits << 1 == roundIncrement)) {
            roundMask |= roundIncrement;
        }

        sig &= ~roundMask;
        zSPtr->signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
        zSPtr->signif = sig;
        return;
    }

    uint32_t sigExtra = extSigPtr[indexWordLo(3)];
    bool doIncrement = UINT32_C(0x80000000) <= sigExtra;

    if (softfloat_round_near_even != softfloat_roundingMode && softfloat_round_near_maxMag != softfloat_roundingMode) {
        doIncrement =
            (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode &&
            sigExtra;
    }

    if (static_cast<uint32_t>(exp - 1) < 0x7FFD) {

        if (sigExtra) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        if (doIncrement) {
            ++sig;

            if (!sig) {
                ++exp;
                sig = static_cast<uint64_t>(INT64_MIN);
            } else {
                sig &= ~static_cast<uint64_t>(!(sigExtra & INT32_MAX) & !!(softfloat_round_near_even == softfloat_roundingMode));
            }
        }

        zSPtr->signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
        zSPtr->signif = sig;
        return;
    }

    if (exp <= 0) {
        bool const isTiny =
            softfloat_tininess_beforeRounding == softfloat_detectTininess ||
            exp < 0 ||
            !doIncrement ||
            sig < UINT64_MAX;
        /** @todo use local variable for output */
        softfloat_shiftRightJam96M(extSigPtr, static_cast<uint8_t>(1 - exp), extSigPtr);
        sigExtra = extSigPtr[indexWordLo(3)];

        if (sigExtra) {
            if (isTiny) {
                softfloat_raiseFlags(softfloat_flag_underflow);
            }

            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        doIncrement =
            softfloat_round_near_even == softfloat_roundingMode || softfloat_round_near_maxMag == softfloat_roundingMode ? static_cast<uint32_t>(INT32_MIN) <= sigExtra :
            (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode && 0 != sigExtra;
        exp = 0;
        sig = static_cast<uint64_t>(extSigPtr[indexWord(3, 2)]) << 32 | extSigPtr[indexWord(3, 1)];

        if (doIncrement) {
            ++sig;
            sig &= ~static_cast<uint64_t>(!(sigExtra & UINT32_C(0x7FFFFFFF)) & !!(softfloat_round_near_even == softfloat_roundingMode));
            exp = 0 != (sig & UINT64_C(0x8000000000000000));
        }

        zSPtr->signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
        zSPtr->signif = sig;
        return;
    }

    if (0x7FFE < exp || (exp == 0x7FFE && sig == UINT64_MAX && doIncrement)) {
        uint64_t const roundMask = 0;
        softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);

        if (
            softfloat_round_near_even == softfloat_roundingMode ||
            softfloat_round_near_maxMag == softfloat_roundingMode ||
            (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode
        ) {
            exp = INT16_MAX;
            sig = static_cast<uint64_t>(INT64_MIN);
        } else {
            exp = INT16_MAX - 1;
            sig = ~roundMask;
        }

        zSPtr->signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
        zSPtr->signif = sig;
        return;
    }
}

}  // namespace softfloat
