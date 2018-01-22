
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

#include "internals.hpp"

#include "specialize.hpp"
#include "softfloat/functions.h"

namespace softfloat {
namespace internals {

float16_t
softfloat_subMagsF16(uint16_t uiA,
                     uint16_t uiB)
{
    uint16_t uiZ;
    int8_t expZ;
    uint16_t sigZ;


    int8_t expA = expF16UI(uiA);
    uint16_t const sigA = fracF16UI(uiA);
    int8_t const expB = expF16UI(uiB);
    uint16_t const sigB = fracF16UI(uiB);
    int8_t expDiff = expA - expB;

    if (!expDiff) {
        if (expA == 0x1F) {
            if (sigA | sigB) {
                return u_as_f_16(softfloat_propagateNaNF16UI(uiA, uiB));
            }

            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_16(defaultNaNF16UI);
        }

        int16_t sigDiff = sigA - sigB;

        if (!sigDiff) {
            return u_as_f_16(packToF16UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
        }

        if (expA) {
            --expA;
        }

        bool signZ = signF16UI(uiA);

        if (sigDiff < 0) {
            signZ = !signZ;
            sigDiff = -sigDiff;
        }

        int8_t shiftDist = softfloat_countLeadingZeros16(static_cast<uint16_t>(sigDiff)) - 5;
        expZ = expA - shiftDist;

        if (expZ < 0) {
            shiftDist = expA;
            expZ = 0;
        }

        return u_as_f_16(packToF16UI(signZ, expZ, static_cast<uint16_t>(sigDiff << shiftDist)));
    }

    {
        bool signZ = signF16UI(uiA);
        uint16_t sigY;
        uint16_t sigX;

        if (expDiff < 0) {
            signZ = !signZ;

            if (expB == 0x1F) {
                return u_as_f_16(sigB ? softfloat_propagateNaNF16UI(uiA, uiB) : packToF16UI(signZ, 0x1F, 0));
            }

            if (expDiff <= -13) {
                uiZ = packToF16UI(signZ, expB, sigB);

                if (0 == (expA | sigA)) {
                    return u_as_f_16(uiZ);
                }

                goto subEpsilon;
            }

            expZ = expA + 19;
            sigX = sigB | 0x0400u;
            sigY = sigA + (expA ? 0x0400u : sigA);
            expDiff = -expDiff;
        } else {
            uiZ = uiA;

            if (expA == 0x1F) {
                return u_as_f_16(sigA ? softfloat_propagateNaNF16UI(uiA, uiB) : uiZ);
            }

            if (!(13 <= expDiff)) {
                expZ = expB + 19;
                sigX = (sigA | 0x0400u);
                sigY = sigB + (expB ? 0x0400u : sigB);
            } else {
                if (!(expB | sigB)) {
                    return u_as_f_16(uiZ);
                }

                goto subEpsilon;
            }
        }

        uint32_t sig32Z = (static_cast<uint32_t>(sigX) << expDiff) - sigY;
        int8_t shiftDist = softfloat_countLeadingZeros32(sig32Z) - 1;
        sig32Z <<= shiftDist;
        expZ -= shiftDist;
        sigZ = sig32Z >> 16;

        if (sig32Z & 0xFFFF) {
            sigZ |= 1;
        } else if (!(sigZ & 0xF) && ((unsigned int)expZ < 0x1E)) {
            sigZ >>= 4;
            return u_as_f_16(packToF16UI(signZ, expZ, sigZ));
        }

        return softfloat_roundPackToF16(signZ, expZ, sigZ);
    }

subEpsilon:
    {
        if (
            softfloat_round_near_even != softfloat_roundingMode &&
            (softfloat_round_minMag == softfloat_roundingMode || (signF16UI(uiZ) ? softfloat_round_max : softfloat_round_min) == softfloat_roundingMode)
        ) {
            --uiZ;
        }

        softfloat_raiseFlags(softfloat_flag_inexact);
        return u_as_f_16(uiZ);
    }
}

}  // namespace internals
}  // namespace softfloat
