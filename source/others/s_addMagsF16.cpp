
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
#include "softfloat/functions.h"

#include "specialize.hpp"

namespace softfloat {

float16_t
softfloat_addMagsF16(uint16_t uiA,
                     uint16_t uiB)
{
    int8_t const expA = expF16UI(uiA);
    uint16_t const sigA = fracF16UI(uiA);
    int8_t const expB = expF16UI(uiB);
    uint16_t const sigB = fracF16UI(uiB);
    int8_t const expDiff = expA - expB;
    uint16_t sigZ;
    int8_t expZ;
    bool signZ;

    if (!expDiff) {
        if (!expA) {
            return u_as_f_16(static_cast<uint16_t>(uiA + sigB));
        }

        if (expA == 0x1F) {
            return u_as_f_16(sigA | sigB ? softfloat_propagateNaNF16UI(uiA, uiB) : uiA);
        }

        signZ = signF16UI(uiA);
        expZ = expA;
        sigZ = 0x0800u + sigA + sigB;

        if (!(sigZ & 1) && (expZ < 0x1E)) {
            return u_as_f_16(packToF16UI(signZ, expZ, static_cast<uint16_t>(sigZ >> 1)));
        }

        sigZ <<= 3;
    } else {
        signZ = signF16UI(uiA);
        int8_t shiftDist;
        uint16_t sigY;
        uint16_t sigX;

        if (expDiff < 0) {
            if (expB == 0x1F) {
                return u_as_f_16(sigB ? softfloat_propagateNaNF16UI(uiA, uiB) : packToF16UI(signZ, 0x1F, 0));
            }

            if (expDiff <= -13) {
                uint16_t uiZ = packToF16UI(signZ, expB, sigB);

                if (!(expA | sigA)) {
                    return u_as_f_16(uiZ);
                }

                softfloat_round_mode const roundingMode = softfloat_roundingMode;

                if (softfloat_round_near_even != roundingMode) {
                    if (roundingMode == (signF16UI(uiZ) ? softfloat_round_min : softfloat_round_max)) {
                        ++uiZ;

                        if (static_cast<uint16_t>(uiZ << 1) == UINT16_C(0xF800)) {
                            softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);
                        }
                    }
                }

                softfloat_raiseFlags(softfloat_flag_inexact);
                return u_as_f_16(uiZ);
            }

            expZ = expB;
            sigX = uint16_t(sigB | 0x0400u);
            sigY = sigA + (expA ? 0x0400u : sigA);
            shiftDist = 19 + expDiff;
        } else {
            uint16_t uiZ = uiA;

            if (expA == 0x1F) {
                return u_as_f_16(sigA ? softfloat_propagateNaNF16UI(uiA, uiB) : uiZ);
            }

            if (13 <= expDiff) {
                if (!(expB | sigB)) {
                    return u_as_f_16(uiZ);
                }

                softfloat_round_mode const roundingMode = softfloat_roundingMode;

                if (roundingMode != softfloat_round_near_even) {
                    if (roundingMode == (signF16UI(uiZ) ? softfloat_round_min : softfloat_round_max)) {
                        ++uiZ;

                        if (static_cast<uint16_t>(uiZ << 1) == UINT16_C(0xF800)) {
                            softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);
                        }
                    }
                }

                softfloat_raiseFlags(softfloat_flag_inexact);
                return u_as_f_16(uiZ);
            } else {
                expZ = expA;
                sigX = sigA | 0x0400u;
                sigY = sigB + (expB ? 0x0400u : sigB);
                shiftDist = 19 - expDiff;
            }
        }

        uint32_t sig32Z = (static_cast<uint32_t>(sigX) << 19) + (static_cast<uint32_t>(sigY) << shiftDist);

        if (sig32Z < 0x40000000) {
            --expZ;
            sig32Z <<= 1;
        }

        sigZ = sig32Z >> 16;

        if (sig32Z & 0xFFFF) {
            sigZ |= 1;
        } else if (!(sigZ & 0xF) && (expZ < 0x1E)) {
            return u_as_f_16(packToF16UI(signZ, expZ, static_cast<uint16_t>(sigZ >> 4)));
        }
    }

    return softfloat_roundPackToF16(signZ, expZ, sigZ);
}

}  // namespace softfloat
