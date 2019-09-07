
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
softfloat_addMagsF16(uint16_t const& uiA,
                     uint16_t const& uiB)
{
    int8_t const expA = get_exp(uiA);
    uint16_t const sigA = get_frac(uiA);
    int8_t const expB = get_exp(uiB);
    uint16_t const sigB = get_frac(uiB);
    int8_t const expDiff = expA - expB;

    if (0 == expDiff) {
        if (0 == expA) {
            return u_as_f(static_cast<uint16_t>(uiA + sigB));
        }

        if (0x1F == expA) {
            return u_as_f(0 != (sigA | sigB) ? propagate_NaN(uiA, uiB) : uiA);
        }

        bool const signZ = is_sign(uiA);
        int8_t const expZ = expA;
        uint16_t const sigZ = 0x0800u + sigA + sigB;

        if (0 == (sigZ & 1) && expZ < 0x1E) {
            return u_as_f(packToF16UI(signZ, expZ, static_cast<uint16_t>(sigZ >> 1)));
        }

        return softfloat_roundPackToF16(signZ, expZ, static_cast<uint16_t>(sigZ << 3));
    }

    bool const signZ = is_sign(uiA);
    int8_t shiftDist;
    uint16_t sigY;
    uint16_t sigX;
    int8_t expZ;
    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    if (expDiff < 0) {
        if (!is_finite(uiB)) {
            return is_NaN(uiB) ? u_as_f(propagate_NaN(uiA, uiB)) : make_signed_inf<float16_t>(signZ);
        }

        if (expDiff <= -13) {
            uint16_t uiZ = packToF16UI(signZ, expB, sigB);

            if (0 != (expA | sigA)) {
                softfloat_raiseFlags(softfloat_flag_inexact);

                if (softfloat_round_near_even != softfloat_roundingMode && (is_sign(uiZ) ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode) {
                    ++uiZ;

                    if (static_cast<uint16_t>(uiZ << 1) == UINT16_C(0xF800)) {
                        softfloat_raiseFlags(softfloat_flag_overflow);
                    }
                }
            }

            return u_as_f(uiZ);
        }

        expZ = expB;
        sigX = static_cast<uint16_t>(sigB | 0x0400u);
        sigY = sigA + (expA ? 0x0400u : sigA);
        shiftDist = 19 + expDiff;
    } else {
        uint16_t uiZ = uiA;

        if (expA == 0x1F) {
            return u_as_f(sigA ? propagate_NaN(uiA, uiB) : uiZ);
        }

        if (13 <= expDiff) {
            if (0 != (expB | sigB)) {
                softfloat_raiseFlags(softfloat_flag_inexact);

                if (softfloat_round_near_even != softfloat_roundingMode && (is_sign(uiZ) ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode) {
                    ++uiZ;

                    if (UINT16_C(0xF800) == static_cast<uint16_t>(uiZ << 1)) {
                        softfloat_raiseFlags(softfloat_flag_overflow);
                    }
                }
            }

            return u_as_f(uiZ);
        }

        expZ = expA;
        sigX = sigA | 0x0400u;
        sigY = sigB + (expB ? 0x0400u : sigB);
        shiftDist = 19 - expDiff;
    }

    uint32_t sig32Z =
        (static_cast<uint32_t>(sigX) << 19) +
        (static_cast<uint32_t>(sigY) << shiftDist);

    if (sig32Z < 0x40000000) {
        --expZ;
        sig32Z <<= 1;
    }

    uint16_t const sigZ = sig32Z >> 16;

    if (0 != (sig32Z & 0xFFFF)) {
        return softfloat_roundPackToF16(signZ, expZ, sigZ | 1u);
    }

    if (0 == (sigZ & 0xF) && expZ < 0x1E) {
        return u_as_f(packToF16UI(signZ, expZ, static_cast<uint16_t>(sigZ >> 4)));
    }

    return softfloat_roundPackToF16(signZ, expZ, sigZ);
}

}  // namespace internals
}  // namespace softfloat
