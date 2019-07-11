
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

namespace {
using namespace softfloat::internals;

static float32_t
addMags(uint32_t const uiA,
        uint32_t const uiB)
{
    assert(!isNaNF32UI(uiA) && !isNaNF32UI(uiB));

    if (isInf32UI(uiA) || isInf32UI(uiB)) {
        /** propagate infinity if operand(s) is infinity */
        return signed_inf_F32(is_sign(uiA));
    }

    static uint32_t const hidden_bit = UINT32_C(1) << 23;
    int16_t const expA = expF32UI(uiA);
    uint32_t const sigA = fracF32UI(uiA);
    int16_t const expB = expF32UI(uiB);
    uint32_t const sigB = fracF32UI(uiB);
    int16_t const expDiff = expA - expB;

    if (0 == expDiff) {
        /* if same exponent, fractional are aligned */
        if (0 == expA) {
            /** if a and b are subnormal(s) or zero(s), result is simple sum, possible carry of high sum bit to exponent is valid */
            return u_as_f_32(uiA + sigB);
        }

        /** add hidden bits to operands if fractional are normal aligned */
        uint32_t const sigZ = (sigA + hidden_bit) + (sigB + hidden_bit);
        bool const signA = is_sign(uiA);

        /** if lowest bit is 0 and exponent allow add carry bit of sum without float class change pack shift right significand with carry bit to exponent */
        if (0 == (sigZ & 1) && expA < 0xFE) {
            return u_as_f_32(packToF32UI(signA, expA, sigZ >> 1));
        }

        return softfloat_roundPackToF32(signA, expA, sigZ << 6);
    }

    /* fractional bits before shift are [22..0] and after shift are [28..6] */
    assert(0 == ((sigA | sigB) & (~UINT32_C(0) << 23)));
    uint32_t sigA_scaled = sigA << 6;
    uint32_t sigB_scaled = sigB << 6;
    /* bit before point is [29] */
    uint32_t const hidden_bit_scaled = hidden_bit << 6;
    assert(0 == ((sigA_scaled | sigB_scaled) & (~UINT32_C(0) << 29)));

    /* unaligned fractional */
    int16_t expZ;

    if (expDiff < 0) {
        /* magnitude b greater than magnitude a */
        expZ = expB;
        /* add hidden bit and shift */
        sigA_scaled = softfloat_shiftRightJam32(sigA_scaled + (expA ? hidden_bit_scaled : sigA_scaled), static_cast<uint16_t>(-expDiff));
        sigB_scaled += hidden_bit_scaled;
    } else {
        /* magnitude a greater than magnitude b */
        expZ = expA;
        /* add hidden bit and shift */
        sigB_scaled = softfloat_shiftRightJam32(sigB_scaled + (expB ? hidden_bit_scaled : sigB_scaled), static_cast<uint16_t>(expDiff));
        sigA_scaled += hidden_bit_scaled;
    }

    /* mantissa bits are [29..0] with 1 digit before point */
    uint32_t const sigZ = sigA_scaled + sigB_scaled;

    /* mantissa bits are [30..0] with up to 2 digits before point */
    /* if high mantissa bit is 0, then shift left to adjust leading mantissa bit to position [30] */
    if (sigZ < 2 * hidden_bit_scaled) {
        return softfloat_roundPackToF32(is_sign(uiA), expZ - 1, sigZ << 1);
    }

    return softfloat_roundPackToF32(is_sign(uiA), expZ, sigZ);
}

static float32_t
subMags(uint32_t const uiA,
        uint32_t const uiB)
{
    assert(!isNaNF32UI(uiA) && !isNaNF32UI(uiB));
    int16_t expA = expF32UI(uiA);
    uint32_t const sigA = fracF32UI(uiA);
    int16_t const expB = expF32UI(uiB);
    uint32_t const sigB = fracF32UI(uiB);

    int16_t const expDiff = expA - expB;
    static int16_t const max_exp = 0xFF;

    if (0 == expDiff) {
        if (max_exp == expA) {
            if (0 != sigA || 0 != sigB) {
                return u_as_f_32(propagate_NaN(uiA, uiB));
            }

            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_32(defaultNaNF32UI);
        }

        int32_t const sigDiff = static_cast<int32_t>(sigA - sigB);

        if (0 == sigDiff) {
            softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
            return signed_zero_F32(softfloat_round_min == softfloat_roundingMode);
        }

        if (expA) {
            --expA;
        }

        bool const signZ = sigDiff < 0 ? !is_sign(uiA) : is_sign(uiA);
        int32_t const sigAbsDiff = sigDiff < 0 ? -sigDiff : sigDiff;
        int8_t const shiftDist = count_leading_zeros(static_cast<uint32_t>(sigAbsDiff)) - 8;
        int16_t const expZ = expA - shiftDist;
        return
            u_as_f_32(expZ < 0 ?
                      packToF32UI(signZ, 0, static_cast<uint32_t>(sigAbsDiff << expA)) :
                      packToF32UI(signZ, expZ, static_cast<uint32_t>(sigAbsDiff << shiftDist)));
    }

    bool const signZ = is_sign(uiA);
    auto const sigA_1 = sigA << 7;
    auto const sigB_1 = sigB << 7;

    if (expDiff < 0) {
        return
            max_exp != expB ?
            softfloat_normRoundPackToF32(!signZ, expB - 1, (sigB_1 | 0x40000000) - softfloat_shiftRightJam32(sigA_1 + (expA ? 0x40000000 : sigA_1), static_cast<uint16_t>(-expDiff))) :
            0 == sigB_1 ?
            signed_inf_F32(!signZ) :
            u_as_f_32(propagate_NaN(uiA, uiB));
    }

    return
        max_exp == expA ?
        u_as_f_32(0 != sigA_1 ? propagate_NaN(uiA, uiB) : uiA) :
        softfloat_normRoundPackToF32(signZ,
                                     expA - 1,
                                     (sigA_1 | 0x40000000) - softfloat_shiftRightJam32(sigB_1 + (expB ? 0x40000000 : sigB_1), static_cast<uint16_t>(expDiff)));
}

}  // namespace

float32_t
f32_add(float32_t const a,
        float32_t const b)
{
    using namespace softfloat::internals;
    uint32_t const uiA = f_as_u_32(a);
    uint32_t const uiB = f_as_u_32(b);

    if (isNaNF32UI(uiA) || isNaNF32UI(uiB)) {
        /** propagate NaN if operand(s) is NaN*/
        return u_as_f_32(propagate_NaN(uiA, uiB));
    }

    return
        (is_sign(uiA) == is_sign(uiB) ? addMags : subMags)(uiA, uiB);
}

float32_t
f32_sub(float32_t const a,
        float32_t const b)
{
    using namespace softfloat::internals;

    uint32_t const uiA = f_as_u_32(a);
    uint32_t const uiB = f_as_u_32(b);

    if (isNaNF32UI(uiA) || isNaNF32UI(uiB)) {
        /** propagate NaN if operand(s) is NaN*/
        return u_as_f_32(propagate_NaN(uiA, uiB));
    }

    return (is_sign(uiA) == is_sign(uiB) ? subMags : addMags)(uiA, uiB);
}
