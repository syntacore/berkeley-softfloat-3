
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
All rights reserved.

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

static inline float32_t
mulAdd(uint32_t const& uiA,
       uint32_t const& uiB,
       uint32_t const& uiC)
{
    if (is_NaN(uiA) || is_NaN(uiB) || is_sNaN(uiC)) {
        return u_as_f_32(propagate_NaN(propagate_NaN(uiA, uiB), uiC));
    }

    bool const signA = is_sign(uiA);
    int16_t expA = get_exp(uiA);
    uint32_t sigA = get_frac(uiA);
    bool const signB = is_sign(uiB);
    int16_t expB = get_exp(uiB);
    uint32_t sigB = get_frac(uiB);
    bool const signC = is_sign(uiC);
    int16_t expC = get_exp(uiC);
    uint32_t sigC = get_frac(uiC);
    bool const signProd = signA != signB;

    static int16_t const max_exp = 0xFF;

    if (max_exp == expA || max_exp == expB) {
        bool const is_product_undefined =
            max_exp == expA ? /* a is inf, check b for zero  */ 0 == expB && 0 == sigB :
            /* expB == max_exp b is inf, check a for zero*/ 0 == expA && 0 == sigA;

        if (is_product_undefined) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_32(propagate_NaN(defaultNaNF32UI, uiC));
        }

        /* product is inf */
        float32_t const uiZ = make_signed_inf<float32_t>(signProd);

        if (expC != max_exp) {
            /* summand c is finite, return product as result */
            return uiZ;
        }

        /* if summand is inf, check for same sign */
        if (is_NaN(uiC) && signProd == signC) {
            /* summands are same sign inf */
            return uiZ;
        }

        /* summands are different sign inf or NaN, undefined sum */
        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_32(propagate_NaN(defaultNaNF32UI, uiC));
    }

    if (is_NaN(uiC)) {
        return u_as_f_32(propagate_NaN(defaultNaNF32UI, uiC));
    }

    if (max_exp == expC) {
        /** if c is infinity while a and b are finite, return c */
        return u_as_f_32(uiC);
    }

    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    if (0 == expA) {
        /* a is zero or subnormal */
        if (0 == sigA) {
            /* a is zero */
            return 0 == (expC | sigC) && signProd != signC ? make_signed_zero<float32_t>(softfloat_round_min == softfloat_roundingMode) : u_as_f_32(uiC);
        }

        exp16_sig32 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (0 == expB) {
        /* b is zero or subnormal */
        if (0 == sigB) {
            /* b is zero */
            return
                0 == (expC | sigC) && signProd != signC ?
                make_signed_zero<float32_t>(softfloat_round_min == softfloat_roundingMode) :
                u_as_f_32(uiC);
        }

        exp16_sig32 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int16_t expProd = expA + expB - 0x7E;
    sigA = (sigA | 0x00800000) << 7;
    sigB = (sigB | 0x00800000) << 7;
    uint64_t sigProd = static_cast<uint64_t>(sigA) * sigB;

    if (sigProd < UINT64_C(0x2000000000000000)) {
        --expProd;
        sigProd <<= 1;
    }

    bool signZ = signProd;

    if (0 == expC) {
        if (0 == sigC) {
            return softfloat_roundPackToF32(signZ, expProd - 1, static_cast<uint32_t>(softfloat_shortShiftRightJam64(sigProd, 31)));
        }

        exp16_sig32 const normExpSig(sigC);
        expC = normExpSig.exp;
        sigC = normExpSig.sig;
    }

    sigC = (sigC | 0x00800000) << 6;

    int16_t const expDiff = expProd - expC;

    if (signProd == signC) {
        int16_t expZ;
        uint32_t sigZ;

        if (expDiff <= 0) {
            expZ = expC;
            sigZ = sigC + static_cast<uint32_t>(softfloat_shiftRightJam64(sigProd, static_cast<uint32_t>(32 - expDiff)));
        } else {
            expZ = expProd;
            uint64_t const sig64Z = sigProd + softfloat_shiftRightJam64(static_cast<uint64_t>(sigC) << 32, static_cast<uint32_t>(expDiff));
            sigZ = static_cast<uint32_t>(softfloat_shortShiftRightJam64(sig64Z, 32));
        }

        if (sigZ < 0x40000000) {
            --expZ;
            sigZ <<= 1;
        }

        return softfloat_roundPackToF32(signZ, expZ, sigZ);
    }

    uint64_t const sig64C = static_cast<uint64_t>(sigC) << 32;
    uint64_t sig64Z;
    int16_t expZ;

    if (expDiff < 0) {
        signZ = signC;
        expZ = expC;
        sig64Z = sig64C - softfloat_shiftRightJam64(sigProd, static_cast<uint32_t>(-expDiff));
    } else if (!expDiff) {
        expZ = expProd;
        sig64Z = sigProd - sig64C;

        if (0 == sig64Z) {
            return make_signed_zero<float32_t>(softfloat_round_min == softfloat_roundingMode);
        }

        if (sig64Z & INT64_MIN) {
            signZ = !signZ;
            sig64Z = static_cast<uint64_t>(-static_cast<int64_t>(sig64Z));
        }
    } else {
        expZ = expProd;
        sig64Z = sigProd - softfloat_shiftRightJam64(sig64C, static_cast<uint32_t>(expDiff));
    }

    int8_t shiftDist = count_leading_zeros(sig64Z) - 1;
    expZ -= shiftDist;
    shiftDist -= 32;

    uint32_t const sigZ =
        shiftDist < 0 ?
        static_cast<uint32_t>(softfloat_shortShiftRightJam64(sig64Z, static_cast<uint8_t>(-shiftDist))) :
        static_cast<uint32_t>(sig64Z) << shiftDist;

    return softfloat_roundPackToF32(signZ, expZ, sigZ);
}

}  // namespace

float32_t
f32_mulAdd(float32_t a,
           float32_t b,
           float32_t c)
{
    return mulAdd(f_as_u(a), f_as_u(b), f_as_u(c));
}

