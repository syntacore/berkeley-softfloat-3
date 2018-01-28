
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

#include "target.hpp"

namespace softfloat {
namespace internals {

extFloat80_t
softfloat_roundPackToExtF80(bool sign,
                            int32_t exp,
                            uint64_t sig,
                            uint64_t sigExtra,
                            uint8_t roundingPrecision)
{
    uint64_t roundIncrement;
    uint64_t roundMask;
    bool doIncrement;
    uint64_extra sig64Extra;
    extFloat80_t uZ;

    bool const roundNearEven = softfloat_round_near_even == softfloat_roundingMode;

    if (roundingPrecision == 80) {
        goto precision80;
    } else if (roundingPrecision == 64) {
        roundIncrement = UINT64_C(0x0000000000000400);
        roundMask = UINT64_C(0x00000000000007FF);
    } else if (roundingPrecision == 32) {
        roundIncrement = UINT64_C(0x0000008000000000);
        roundMask = UINT64_C(0x000000FFFFFFFFFF);
    } else {
        goto precision80;
    }

    sig |= (sigExtra != 0);

    if (!roundNearEven && (softfloat_roundingMode != softfloat_round_near_maxMag)) {
        roundIncrement = (softfloat_roundingMode == (sign ? softfloat_round_min : softfloat_round_max)) ? roundMask : 0;
    }

    uint64_t roundBits = sig & roundMask;

    if (0x7FFD <= static_cast<uint32_t>(exp - 1)) {
        if (exp <= 0) {
            bool isTiny = softfloat_tininess_beforeRounding == softfloat_detectTininess || exp < 0 || sig <= static_cast<uint64_t>(sig + roundIncrement);
            sig = softfloat_shiftRightJam64(sig, 1u - exp);
            roundBits = sig & roundMask;

            if (isTiny && roundBits) {
                softfloat_raiseFlags(softfloat_flag_underflow);
            }

            if (roundBits) {
                softfloat_raiseFlags(softfloat_flag_inexact);
            }

            sig += roundIncrement;
            exp = ((sig & UINT64_C(0x8000000000000000)) != 0);
            roundIncrement = roundMask + 1;

            if (roundNearEven && (roundBits << 1 == roundIncrement)) {
                roundMask |= roundIncrement;
            }

            sig &= ~roundMask;
            uZ.signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
            uZ.signif = sig;
            return uZ;
        }

        if (0x7FFE < exp || (0x7FFE == exp  && static_cast<uint64_t>(sig + roundIncrement) < sig) ) {
            goto overflow;
        }
    }

    if (roundBits) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    sig = static_cast<uint64_t>(sig + roundIncrement);

    if (sig < roundIncrement) {
        ++exp;
        sig = UINT64_C(0x8000000000000000);
    }

    roundIncrement = roundMask + 1;

    if (roundNearEven && (roundBits << 1 == roundIncrement)) {
        roundMask |= roundIncrement;
    }

    sig &= ~roundMask;

    if (!sig) {
        exp = 0;
    }

    uZ.signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
    uZ.signif = sig;
    return uZ;

precision80:
    doIncrement = UINT64_C(0x8000000000000000) <= sigExtra;

    if (!roundNearEven && (softfloat_roundingMode != softfloat_round_near_maxMag)) {
        doIncrement = (softfloat_roundingMode == (sign ? softfloat_round_min : softfloat_round_max)) && sigExtra;
    }

    if (0x7FFD <= static_cast<uint32_t>(exp - 1)) {
        if (exp <= 0) {
            bool isTiny = softfloat_detectTininess == softfloat_tininess_beforeRounding || exp < 0 || !doIncrement || sig < UINT64_C(0xFFFFFFFFFFFFFFFF);
            sig64Extra = softfloat_shiftRightJam64Extra(sig, sigExtra, 1u - exp);
            sig = sig64Extra.v;
            sigExtra = sig64Extra.extra;

            if (isTiny && sigExtra) {
                softfloat_raiseFlags(softfloat_flag_underflow);
            }

            if (sigExtra) {
                softfloat_raiseFlags(softfloat_flag_inexact);
            }

            doIncrement = UINT64_C(0x8000000000000000) <= sigExtra;

            if (!roundNearEven && (softfloat_roundingMode != softfloat_round_near_maxMag)) {
                doIncrement = (softfloat_roundingMode == (sign ? softfloat_round_min : softfloat_round_max)) && sigExtra;
            }

            exp = 0;

            if (doIncrement) {
                ++sig;
                sig &= ~static_cast<uint64_t>(!(sigExtra & UINT64_C(0x7FFFFFFFFFFFFFFF)) & roundNearEven);
                exp = ((sig & UINT64_C(0x8000000000000000)) != 0);
            }

            uZ.signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
            uZ.signif = sig;
            return uZ;
        }

        if (0x7FFE < exp || (exp == 0x7FFE && sig == UINT64_C(0xFFFFFFFFFFFFFFFF) && doIncrement)) {
            roundMask = 0;
overflow:
            softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);

            if (roundNearEven || softfloat_roundingMode == softfloat_round_near_maxMag || softfloat_roundingMode == (sign ? softfloat_round_min : softfloat_round_max)) {
                exp = 0x7FFF;
                sig = UINT64_C(0x8000000000000000);
            } else {
                exp = 0x7FFE;
                sig = ~roundMask;
            }

            uZ.signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
            uZ.signif = sig;
            return uZ;
        }
    }

    if (sigExtra) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    if (doIncrement) {
        ++sig;

        if (!sig) {
            ++exp;
            sig = UINT64_C(0x8000000000000000);
        } else {
            sig &= ~static_cast<uint64_t>(!(sigExtra & UINT64_C(0x7FFFFFFFFFFFFFFF)) & roundNearEven);
        }
    } else {
        if (!sig) {
            exp = 0;
        }
    }

    uZ.signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp));
    uZ.signif = sig;
    return uZ;
}

}  // namespace internals
}  // namespace softfloat
