
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

float32_t
softfloat_subMagsF32(uint32_t uiA,
                     uint32_t uiB)
{
    int16_t expA = expF32UI(uiA);
    uint32_t sigA = fracF32UI(uiA);
    int16_t const expB = expF32UI(uiB);
    uint32_t sigB = fracF32UI(uiB);

    int16_t expDiff = expA - expB;

    if (0 == expDiff) {
        if (0xFF == expA) {
            if (0 != sigA || 0 != sigB) {
                return u_as_f_32(softfloat_propagateNaNF32UI(uiA, uiB));
            }

            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_32(defaultNaNF32UI);
        }

        int32_t sigDiff = static_cast<int32_t>(sigA - sigB);

        if (0 == sigDiff) {
            return signed_zero_F32(softfloat_round_min == softfloat_roundingMode);
        }

        if (expA) {
            --expA;
        }

        bool const signZ = sigDiff < 0 ? !signF32UI(uiA) : signF32UI(uiA);
        int32_t const sigAbsDiff = sigDiff < 0 ? -sigDiff : sigDiff;
        int8_t const shiftDist = softfloat_countLeadingZeros32(static_cast<uint32_t>(sigAbsDiff)) - 8;
        int16_t const expZ = expA - shiftDist;
        return
            u_as_f_32(
                expZ < 0 ? packToF32UI(signZ, 0, static_cast<uint32_t>(sigAbsDiff << expA)) :
                packToF32UI(signZ, expZ, static_cast<uint32_t>(sigAbsDiff << shiftDist)));
    }

    bool const signZ = signF32UI(uiA);
    sigA <<= 7;
    sigB <<= 7;

    if (expDiff < 0) {
        return
            0xFF != expB ? softfloat_normRoundPackToF32(!signZ,
                    expB - 1,
                    (sigB | 0x40000000) -
                    softfloat_shiftRightJam32(sigA + (expA ? 0x40000000 : sigA), static_cast<uint16_t>(-expDiff))) :
            0 == sigB ? signed_inf_F32(!signZ) :
            u_as_f_32(softfloat_propagateNaNF32UI(uiA, uiB));
    }

    return
        0xFF == expA ? u_as_f_32(0 != sigA ? softfloat_propagateNaNF32UI(uiA, uiB) : uiA) :
        softfloat_normRoundPackToF32(signZ,
                                     expA - 1,
                                     (sigA | 0x40000000) -
                                     softfloat_shiftRightJam32(sigB + (expB ? 0x40000000 : sigB), static_cast<uint16_t>(expDiff)));
}

}  // namespace internals
}  // namespace softfloat
