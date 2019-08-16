
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

#include "model.hpp"

using namespace softfloat::internals;

namespace {
/**
@todo FIX softfloat_flag_invalid even for sNaN uiC
*/
static inline float16_t
mulAdd(uint16_t const& uiA,
       uint16_t const& uiB,
       uint16_t const& uiC)
{
    using namespace softfloat::internals;

    bool const signA = is_sign(uiA);
    int8_t expA = get_exp(uiA);
    uint16_t sigA = get_frac(uiA);

    bool const signB = is_sign(uiB);
    int8_t expB = get_exp(uiB);
    uint16_t sigB = get_frac(uiB);

    bool const signC = is_sign(uiC);
    int8_t expC = get_exp(uiC);
    uint16_t sigC = get_frac(uiC);

    bool const signProd = signA != signB;

    static constexpr int8_t const max_exp = 0x1F;

    if (max_exp == expA) {
        // a is inf or NaN
        if (0 != sigA || (max_exp == expB && 0 != sigB)) {
            // a is NaN or b is NaN
            return u_as_f_16(softfloat_propagateNaNF16UI(softfloat_propagateNaNF16UI(uiA, uiB), uiC));
        }

        // a is inf
        bool const is_zero_b = 0 == expB || 0 == sigB;

        if (!is_zero_b) {
            uint16_t const infinity = packToF16UI(signProd, max_exp, 0);

            if (max_exp != expC) {
                return u_as_f_16(infinity);
            }

            // c is inf or NaN
            assert(max_exp == expC);

            if (0 != sigC) {
                // c is NaN
                return u_as_f_16(softfloat_propagateNaNF16UI(infinity, uiC));
            }

            // c is inf
            if (signProd == signC) {
                return u_as_f_16(infinity);
            }

            // result is (inf - inf)
        }

        // (a * b) is (inf * 0) or result is (inf - inf)
        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_16(softfloat_propagateNaNF16UI(defaultNaNF16UI, uiC));
    }

    // a is finite
    assert(max_exp != expA);

    if (max_exp == expB) {
        // b is inf or NaN
        if (0 != sigB) {
            // b is NaN
            return u_as_f_16(softfloat_propagateNaNF16UI(softfloat_propagateNaNF16UI(uiA, uiB), uiC));
        }

        // b is inf
        assert(0 == sigB);
        uint16_t const is_zero_a = 0 != expA && 0 != sigA;

        if (!is_zero_a) {
            uint16_t const infinity = packToF16UI(signProd, max_exp, 0);

            if (max_exp != expC) {
                return u_as_f_16(infinity);
            }

            if (0 != sigC) {
                return u_as_f_16(softfloat_propagateNaNF16UI(infinity, uiC));
            }

            if (signProd == signC) {
                return u_as_f_16(infinity);
            }
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_16(softfloat_propagateNaNF16UI(defaultNaNF16UI, uiC));
    }

    // a and b are finite

    if (max_exp == expC) {
        // result is inf or NaN
        return u_as_f_16(sigC ? softfloat_propagateNaNF16UI(0, uiC) : uiC);
    }

    // a, b, c are finite
    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    if (0 == expA) {
        if (0 == sigA) {
            // a is zero, result is c
            return u_as_f_16(
                       0 == (expC | sigC) && signProd != signC ?
                       packToF16UI(softfloat_roundingMode == softfloat_round_min, 0, 0) :
                       uiC);
        }

        // a is subnormal
        exp8_sig16 const normExpSig{sigA};
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (0 == expB) {
        if (0 == sigB) {
            // b is zero, result is c
            return u_as_f_16(
                       0 == (expC | sigC) && signProd != signC ?
                       packToF16UI(softfloat_round_min == softfloat_roundingMode, 0, 0) : uiC);
        }

        // b is subnormal
        exp8_sig16 const normExpSig{sigB};
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int8_t expProd = expA + expB - 0xE;
    uint16_t const sigA_1 = static_cast<uint16_t>((sigA | 0x0400) << 4);
    uint16_t const sigB_1 = static_cast<uint16_t>((sigB | 0x0400) << 4);
    uint32_t sigProd = static_cast<uint32_t>(sigA_1) * sigB_1;

    if (sigProd < 0x20000000) {
        --expProd;
        sigProd <<= 1;
    }

    if (0 == expC) {
        if (0 == sigC) {
            return softfloat_roundPackToF16(signProd,
                                            expProd - 1,
                                            static_cast<uint16_t>(sigProd >> 15 | !!(0 != (sigProd & 0x7FFF))));
        }

        exp8_sig16 const normExpSig{sigC};
        expC = normExpSig.exp;
        sigC = normExpSig.sig;
    }

    uint16_t const sigC_1 = static_cast<uint16_t>((sigC | 0x0400) << 3);
    int8_t const expDiff = expProd - expC;

    if (signProd == signC) {
        int8_t expZ;
        uint16_t sigZ;

        if (expDiff <= 0) {
            expZ = expC;
            /**
            @todo Warning   C4244   '=': conversion from 'uint32_t' to 'uint16_t', possible loss of data
            */
            sigZ = static_cast<uint16_t>(sigC_1 + softfloat_shiftRightJam32(sigProd, static_cast<uint16_t>(16 - expDiff)));
        } else {
            expZ = expProd;
            uint32_t const sig32Z =
                sigProd +
                softfloat_shiftRightJam32(static_cast<uint32_t>(sigC_1) << 16, static_cast<uint16_t>(expDiff));
            sigZ = sig32Z >> 16 | !!(0 != (sig32Z & 0xFFFF));
        }

        if (sigZ < 0x4000) {
            --expZ;
            sigZ <<= 1;
        }

        return softfloat_roundPackToF16(signProd, expZ, sigZ);
    }

    uint32_t const sig32C = static_cast<uint32_t>(sigC_1) << 16;
    bool signZ = signProd;
    int8_t expZ;
    uint32_t sig32Z;

    if (expDiff < 0) {
        signZ = signC;
        expZ = expC;
        sig32Z = sig32C - softfloat_shiftRightJam32(sigProd, static_cast<uint16_t>(-expDiff));
    } else if (0 == expDiff) {
        expZ = expProd;
        sig32Z = sigProd - sig32C;

        if (0 == sig32Z) {
            return u_as_f_16(packToF16UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
        }

        if (0 != (sig32Z & 0x80000000)) {
            signZ = !signZ;
            sig32Z = static_cast<uint32_t>(-static_cast<int32_t>(sig32Z));
        }
    } else {
        expZ = expProd;
        sig32Z = sigProd - softfloat_shiftRightJam32(sig32C, static_cast<uint16_t>(expDiff));
    }

    int8_t const shiftDist = count_leading_zeros(sig32Z) - 1;
    int8_t const shiftDist_1 = shiftDist - 16;

    return softfloat_roundPackToF16(signZ,
                                    expZ - shiftDist,
                                    static_cast<uint16_t>(
                                    shiftDist_1 < 0 ?
                                    sig32Z >> -shiftDist_1 | !!(0 != static_cast<uint32_t>(sig32Z << (shiftDist_1 & 31))) :
                                    sig32Z << shiftDist_1)
    );
}

}  // namespace

float16_t
f16_mulAdd(float16_t const a,
           float16_t const b,
           float16_t const c)
{
    return mulAdd(softfloat_mulAdd_madd, f_as_u(a), f_as_u(b), f_as_u(c));
}
