
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

#include "internals.h"

#include "specialize.h"
#include "softfloat/functions.h"

float32_t softfloat_subMagsF32(uint32_t uiA, uint32_t uiB)
{
    int32_t sigDiff;
    bool signZ;
    int16_t expZ;

    int16_t expA = expF32UI(uiA);
    uint32_t sigA = fracF32UI(uiA);
    int16_t const expB = expF32UI(uiB);
    uint32_t sigB = fracF32UI(uiB);

    int16_t expDiff = expA - expB;
    if (!expDiff) {
        int8_t shiftDist;
        if (expA == 0xFF) {
            if (sigA | sigB) {
                return u_as_f_32(softfloat_propagateNaNF32UI(uiA, uiB));
            } else {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return u_as_f_32(defaultNaNF32UI);
            }
        }
        sigDiff = sigA - sigB;
        if (!sigDiff) {
            return signed_zero_F32(softfloat_round_min == softfloat_roundingMode);
        }
        if (expA) {
            --expA;
        }
        signZ = signF32UI(uiA);
        if (sigDiff < 0) {
            signZ = !signZ;
            sigDiff = -sigDiff;
        }
        shiftDist = softfloat_countLeadingZeros32(sigDiff) - 8;
        expZ = expA - shiftDist;
        if (expZ < 0) {
            shiftDist = expA;
            expZ = 0;
        }
        return u_as_f_32(packToF32UI(signZ, expZ, sigDiff << shiftDist));
    } else {
        uint32_t sigX, sigY;
        signZ = signF32UI(uiA);
        sigA <<= 7;
        sigB <<= 7;
        if (expDiff < 0) {
            signZ = !signZ;
            if (expB == 0xFF) {
                if (sigB) {
                    return u_as_f_32(softfloat_propagateNaNF32UI(uiA, uiB));
                }
                return signed_inf_F32(signZ);
            }
            expZ = expB - 1;
            sigX = sigB | 0x40000000;
            sigY = sigA + (expA ? 0x40000000 : sigA);
            expDiff = -expDiff;
        } else {
            if (expA == 0xFF) {
                if (sigA) {
                    return u_as_f_32(softfloat_propagateNaNF32UI(uiA, uiB));
                }
                return u_as_f_32(uiA);
            }
            expZ = expA - 1;
            sigX = sigA | 0x40000000;
            sigY = sigB + (expB ? 0x40000000 : sigB);
        }
        return
            softfloat_normRoundPackToF32(signZ, expZ, sigX - softfloat_shiftRightJam32(sigY, expDiff));
    }
}
