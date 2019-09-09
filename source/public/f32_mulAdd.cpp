
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

float32_t
f32_mulAdd(float32_t const a,
           float32_t const b,
           float32_t const c)
{
    using namespace softfloat::internals;

    if (is_NaN(a) || is_NaN(b) || is_sNaN(c)) {
        return propagate_NaN(propagate_NaN(a, b), c);
    }

    bool const signA = is_sign(a);
    bool const signB = is_sign(b);
    bool const signC = is_sign(c);
    bool const signProd = signA != signB;

    if (is_inf(a) || is_inf(b)) {
        /* a or b is inf, product is inf or undefined, check other operand for zero */
        bool const is_product_undefined = is_inf(a) ? is_zero(b) : is_zero(a);

        if (is_product_undefined) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return propagate_NaN(u_as_f(defaultNaNF32UI), c);
        }

        if (is_NaN(c)) {
            return propagate_NaN(u_as_f(defaultNaNF32UI), c);
        }

        if (is_finite(c) || signProd == signC) {
            /* if summand c is finite or same sign return product as result*/
            return make_signed_inf<float32_t>(signProd);
        }

        /* summands are different sign inf, undefined sum */
        softfloat_raiseFlags(softfloat_flag_invalid);
        return propagate_NaN(u_as_f(defaultNaNF32UI), c);
    }

    if (is_NaN(c)) {
        return propagate_NaN(u_as_f(defaultNaNF32UI), c);
    }

    if (is_inf(c)) {
        /** if c is infinity while a and b are finite, return c */
        return c;
    }

    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    if (is_zero(a) || is_zero(b)) {
        return
            is_zero(c) && signProd != signC ?
            make_signed_zero<float32_t>(softfloat_round_min == softfloat_roundingMode) :
            c;
    }

    int16_t expA = get_exp(a);
    uint32_t sigA = get_frac(a);

    if (0 == expA) {
        exp16_sig32 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int16_t expB = get_exp(b);
    uint32_t sigB = get_frac(b);

    if (0 == expB) {
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

    if (is_zero(c)) {
        return round_pack_to_F32(signZ, expProd - 1, static_cast<uint32_t>(short_shift_right_jam_64(sigProd, 31)));
    }

    int16_t expC = get_exp(c);
    uint32_t sigC = get_frac(c);

    if (0 == expC) {
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
            sigZ = sigC + static_cast<uint32_t>(shift_right_jam_64(sigProd, static_cast<uint32_t>(32 - expDiff)));
        } else {
            expZ = expProd;
            uint64_t const sig64Z = sigProd + shift_right_jam_64(static_cast<uint64_t>(sigC) << 32, static_cast<uint32_t>(expDiff));
            sigZ = static_cast<uint32_t>(short_shift_right_jam_64(sig64Z, 32));
        }

        if (sigZ < 0x40000000) {
            --expZ;
            sigZ <<= 1;
        }

        return round_pack_to_F32(signZ, expZ, sigZ);
    }

    uint64_t const sig64C = static_cast<uint64_t>(sigC) << 32;
    uint64_t sig64Z;
    int16_t expZ;

    if (expDiff < 0) {
        signZ = signC;
        expZ = expC;
        sig64Z = sig64C - shift_right_jam_64(sigProd, static_cast<uint32_t>(-expDiff));
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
        sig64Z = sigProd - shift_right_jam_64(sig64C, static_cast<uint32_t>(expDiff));
    }

    int8_t shiftDist = count_leading_zeros(sig64Z) - 1;
    expZ -= shiftDist;
    shiftDist -= 32;

    uint32_t const sigZ =
        shiftDist < 0 ?
        static_cast<uint32_t>(short_shift_right_jam_64(sig64Z, static_cast<uint8_t>(-shiftDist))) :
        static_cast<uint32_t>(sig64Z) << shiftDist;

    return round_pack_to_F32(signZ, expZ, sigZ);
}

