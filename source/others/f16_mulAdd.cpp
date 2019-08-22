
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

/**
@todo FIX softfloat_flag_invalid even for sNaN uiC
*/
float16_t
f16_mulAdd(float16_t const a,
           float16_t const b,
           float16_t const c)
{
    using namespace softfloat::internals;

    uint16_t const uiA = f_as_u(a);
    uint16_t const uiB = f_as_u(b);
    uint16_t const uiC = f_as_u(c);

    if (is_NaN(uiA) || is_NaN(uiB) || is_sNaN(uiC)) {
        return to_float(propagate_NaN(propagate_NaN(uiA, uiB), uiC));
    }

    bool const signA = is_sign(uiA);
    bool const signB = is_sign(uiB);
    bool const signC = is_sign(uiC);
    bool const signProd = signA != signB;

    if (is_inf(uiA) || is_inf(uiB)) {
        /* a or b is inf, product is inf or undefined, check other operand for zero */
        bool const is_product_undefined = is_inf(uiA) ? is_zero(uiB) : is_zero(uiA);

        if (is_product_undefined) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return to_float(propagate_NaN(defaultNaNF16UI, uiC));
        }

        if (is_NaN(uiC)) {
            return to_float(propagate_NaN(defaultNaNF16UI, uiC));
        }

        if (is_finite(uiC) || signProd == signC) {
            /* if summand c is finite or same sign return product as result*/
            return make_signed_inf<float16_t>(signProd);
        }

        /* summands are different sign inf, undefined sum */
        softfloat_raiseFlags(softfloat_flag_invalid);
        return to_float(propagate_NaN(defaultNaNF16UI, uiC));
    }

    if (is_NaN(uiC)) {
        return to_float(propagate_NaN(defaultNaNF16UI, uiC));
    }

    if (is_inf(uiC)) {
        /** if c is infinity while a and b are finite, return c */
        return to_float(uiC);
    }

    // a, b, c are finite
    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
    int8_t expA = get_exp(uiA);
    uint16_t sigA = get_frac(uiA);

    if (0 == expA) {
        /* a is zero or subnormal */
        if (0 == sigA) {
            // a is zero, result is c
            return is_zero(uiC) && signProd != signC ?
                make_signed_zero<float16_t>(softfloat_round_min == softfloat_roundingMode) :
                to_float(uiC);
        }

        // a is subnormal
        exp8_sig16 const normExpSig{sigA};
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int8_t expB = get_exp(uiB);
    uint16_t sigB = get_frac(uiB);

    if (0 == expB) {
        if (0 == sigB) {
            // b is zero, result is c
            return is_zero(uiC) && signProd != signC ?
                make_signed_zero<float16_t>(softfloat_round_min == softfloat_roundingMode) :
                to_float(uiC);
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

    int8_t expC = get_exp(uiC);
    uint16_t sigC = get_frac(uiC);

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
            return to_float(packToF16UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
        }

        if (0 != (sig32Z & 0x80000000)) {
            signZ = !signZ;
            sig32Z = static_cast<uint32_t>(-static_cast<int32_t>(sig32Z));
        }
    } else /*if (0 < expDiff)*/ {
        expZ = expProd;
        sig32Z = sigProd - softfloat_shiftRightJam32(sig32C, static_cast<uint16_t>(expDiff));
    }

    int8_t const shiftDist = count_leading_zeros(sig32Z) - 1;
    int8_t const shiftDist_1 = shiftDist - 16;
    uint16_t const sigZ =
        shiftDist_1 < 0 ?
        sig32Z >> -shiftDist_1 | !!(0 != static_cast<uint32_t>(sig32Z << (shiftDist_1 & 31))) :
        sig32Z << shiftDist_1;

    return softfloat_roundPackToF16(signZ, expZ - shiftDist, sigZ);
}
