
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

#include "model.hpp"

namespace softfloat {
namespace internals {

float16_t
sub_magnitudes(uint16_t const& uiA,
               uint16_t const& uiB)
{
    using namespace softfloat::internals;

    int8_t expA = get_exp(uiA);
    uint16_t const sigA = get_frac(uiA);
    int8_t const expB = get_exp(uiB);
    uint16_t const sigB = get_frac(uiB);
    int8_t expDiff = expA - expB;
    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    if (0 == expDiff) {
        if (expA == 0x1F) {
            if (0 != (sigA | sigB)) {
                return u_as_f(propagate_NaN(uiA, uiB));
            }

            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f(defaultNaNF16UI);
        }

        int16_t sigDiff = sigA - sigB;

        if (0 == sigDiff) {
            return u_as_f(pack_to_F16_UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
        }

        if (0 != expA) {
            --expA;
        }

        bool signZ = is_sign(uiA);

        if (sigDiff < 0) {
            signZ = !signZ;
            sigDiff = -sigDiff;
        }

        int8_t const shiftDist = count_leading_zeros(static_cast<uint16_t>(sigDiff)) - 5;
        int8_t const expZ = expA - shiftDist;

        if (expZ < 0) {
            return u_as_f(pack_to_F16_UI(signZ, 0, static_cast<uint16_t>(sigDiff << expA)));
        }

        return u_as_f(pack_to_F16_UI(signZ, expZ, static_cast<uint16_t>(sigDiff << shiftDist)));
    }

    bool signZ = is_sign(uiA);
    uint16_t sigY;
    uint16_t sigX;
    int8_t expZ;

    if (expDiff < 0) {
        signZ = !signZ;

        if (expB == 0x1F) {
            return u_as_f(sigB ? propagate_NaN(uiA, uiB) : pack_to_F16_UI(signZ, 0x1F, 0));
        }

        if (expDiff <= -13) {
            uint16_t const uiZ_1 = pack_to_F16_UI(signZ, expB, sigB);

            if (0 != (expA | sigA)) {
                softfloat_raiseFlags(softfloat_flag_inexact);

                if (
                    softfloat_round_near_even != softfloat_roundingMode &&
                    (softfloat_round_minMag == softfloat_roundingMode || (is_sign(uiZ_1) ? softfloat_round_max : softfloat_round_min) == softfloat_roundingMode)
                ) {
                    return u_as_f(uint16_t(uiZ_1 - 1u));
                }
            }

            return u_as_f(uiZ_1);
        }

        expZ = expA + 19;
        sigX = sigB | 0x0400u;
        sigY = sigA + (expA ? 0x0400u : sigA);
        expDiff = -expDiff;
    } else {
        if (expA == 0x1F) {
            return u_as_f(sigA ? propagate_NaN(uiA, uiB) : uiA);
        }

        if (13 <= expDiff) {
            if (0 != (expB | sigB)) {
                softfloat_raiseFlags(softfloat_flag_inexact);

                if (
                    softfloat_round_near_even != softfloat_roundingMode &&
                    (softfloat_round_minMag == softfloat_roundingMode || (is_sign(uiA) ? softfloat_round_max : softfloat_round_min) == softfloat_roundingMode)
                ) {
                    return u_as_f(uint16_t(uiA - 1u));
                }
            }

            return u_as_f(uiA);
        }

        expZ = expB + 19;
        sigX = sigA | 0x0400u;
        sigY = sigB + (expB ? 0x0400u : sigB);
    }

    uint32_t const sig32Z = (static_cast<uint32_t>(sigX) << expDiff) - sigY;
    int8_t const shiftDist = count_leading_zeros(sig32Z) - 1;
    uint32_t const sig32Z_1 = sig32Z << shiftDist;
    int8_t const expZ_1 = expZ - shiftDist;
    uint16_t const sigZ = sig32Z_1 >> 16;

    if (0 != (sig32Z_1 & 0xFFFFu)) {
        return round_pack_to_F16(signZ, expZ_1, sigZ | 1u);
    }

    if (0 == (sigZ & 0xFu) && static_cast<unsigned>(expZ_1) < 0x1Eu) {
        return u_as_f(pack_to_F16_UI(signZ, expZ_1, static_cast<uint16_t>(sigZ >> 4)));
    }

    return round_pack_to_F16(signZ, expZ_1, sigZ);
}

}  // namespace internals
}  // namespace softfloat
