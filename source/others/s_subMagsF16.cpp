
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

#include "target.hpp"

namespace softfloat {
namespace internals {

float16_t
softfloat_subMagsF16(uint16_t const uiA,
                     uint16_t const uiB)
{
    int8_t expA = expF16UI(uiA);
    uint16_t const sigA = fracF16UI(uiA);
    int8_t const expB = expF16UI(uiB);
    uint16_t const sigB = fracF16UI(uiB);
    int8_t expDiff = expA - expB;

    if (0 == expDiff) {
        if (expA == 0x1F) {
            if (0 != (sigA | sigB)) {
                return u_as_f_16(softfloat_propagateNaNF16UI(uiA, uiB));
            } else {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return u_as_f_16(defaultNaNF16UI);
            }
        } else {
            int16_t sigDiff = sigA - sigB;

            if (0 == sigDiff) {
                return u_as_f_16(packToF16UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
            } else {
                if (0 != expA) {
                    --expA;
                }

                bool signZ = signF16UI(uiA);

                if (sigDiff < 0) {
                    signZ = !signZ;
                    sigDiff = -sigDiff;
                }

                int8_t shiftDist = softfloat_countLeadingZeros16(static_cast<uint16_t>(sigDiff)) - 5;
                int8_t const expZ = expA - shiftDist;

                if (expZ < 0) {
                    return u_as_f_16(packToF16UI(signZ, 0, static_cast<uint16_t>(sigDiff << expA)));
                }

                return u_as_f_16(packToF16UI(signZ, expZ, static_cast<uint16_t>(sigDiff << shiftDist)));
            }
        }
    } else {
        bool signZ = signF16UI(uiA);
        uint16_t sigY;
        uint16_t sigX;
        int8_t expZ;

        if (expDiff < 0) {
            signZ = !signZ;

            if (expB == 0x1F) {
                return u_as_f_16(sigB ? softfloat_propagateNaNF16UI(uiA, uiB) : packToF16UI(signZ, 0x1F, 0));
            } else if (expDiff <= -13) {
                uint16_t const uiZ_1 = packToF16UI(signZ, expB, sigB);

                if (0 != (expA | sigA)) {
                    softfloat_raiseFlags(softfloat_flag_inexact);

                    if (
                        softfloat_round_near_even != softfloat_roundingMode &&
                        (softfloat_round_minMag == softfloat_roundingMode || (signF16UI(uiZ_1) ? softfloat_round_max : softfloat_round_min) == softfloat_roundingMode)
                    ) {
                        return u_as_f_16(uiZ_1 - 1u);
                    }
                }

                return u_as_f_16(uiZ_1);
            } else {
                expZ = expA + 19;
                sigX = sigB | 0x0400u;
                sigY = sigA + (expA ? 0x0400u : sigA);
                expDiff = -expDiff;
            }
        } else {
            if (expA == 0x1F) {
                return u_as_f_16(sigA ? softfloat_propagateNaNF16UI(uiA, uiB) : uiA);
            } else if (13 <= expDiff) {
                if (0 != (expB | sigB)) {
                    softfloat_raiseFlags(softfloat_flag_inexact);

                    if (
                        softfloat_round_near_even != softfloat_roundingMode &&
                        (softfloat_round_minMag == softfloat_roundingMode || (signF16UI(uiA) ? softfloat_round_max : softfloat_round_min) == softfloat_roundingMode)
                    ) {
                        return u_as_f_16(uiA - 1u);
                    }
                }

                return u_as_f_16(uiA);
            } else {
                expZ = expB + 19;
                sigX = sigA | 0x0400u;
                sigY = sigB + (expB ? 0x0400u : sigB);
            }
        }

        uint32_t const sig32Z = (static_cast<uint32_t>(sigX) << expDiff) - sigY;
        int8_t const shiftDist = softfloat_countLeadingZeros32(sig32Z) - 1;
        uint32_t const sig32Z_1 = sig32Z << shiftDist;
        int8_t const expZ_1 = expZ - shiftDist;
        uint16_t const sigZ = sig32Z_1 >> 16;

        if (0 != (sig32Z_1 & 0xFFFFu)) {
            return softfloat_roundPackToF16(signZ, expZ_1, sigZ | 1u);
        } else if (0 == (sigZ & 0xFu) && static_cast<unsigned>(expZ_1) < 0x1Eu) {
            return u_as_f_16(packToF16UI(signZ, expZ_1, static_cast<uint16_t>(sigZ >> 4)));
        } else {
            return softfloat_roundPackToF16(signZ, expZ_1, sigZ);
        }
    }
}

}  // namespace internals
}  // namespace softfloat
