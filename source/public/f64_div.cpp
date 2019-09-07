
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

float64_t
f64_div(float64_t const a,
        float64_t const b)
{
    using namespace softfloat::internals;
    bool const signA = is_sign(a);
    bool const signB = is_sign(b);
    bool const signZ = signA != signB;

    if (is_NaN(a) || is_NaN(b)) {
        return propagate_NaN(a, b);
    }

    if (is_inf(a)) {
        if (is_inf(b)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f(defaultNaNF64UI);
        }

        return make_signed_inf<float64_t>(signZ);
    }

    if (is_inf(b)) {
        return make_signed_zero<float64_t>(signZ);
    }

    if (is_zero(b)) {
        if (is_zero(a)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f(defaultNaNF64UI);
        }

        softfloat_raiseFlags(softfloat_flag_infinite);
        return make_signed_inf<float64_t>(signZ);
    }

    if (is_zero(a)) {
        return make_signed_zero<float64_t>(signZ);
    }

    int16_t expB = get_exp(b);
    uint64_t sigB = get_frac(b);

    if (0 == expB) {
        exp16_sig64 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int16_t expA = get_exp(a);
    uint64_t sigA = get_frac(a);

    if (0 == expA) {
        exp16_sig64 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int16_t expZ = expA - expB + 0x3FE;
    sigA |= UINT64_C(0x0010000000000000);
    sigB |= UINT64_C(0x0010000000000000);

    if (sigA < sigB) {
        --expZ;
        sigA <<= 11;
    } else {
        sigA <<= 10;
    }

    sigB <<= 11;
    uint32_t const recip32 = softfloat_approxRecip32_1(sigB >> 32) - 2;
    uint32_t const sig32Z = (static_cast<uint32_t>(sigA >> 32) * static_cast<uint64_t>(recip32)) >> 32;
    uint32_t doubleTerm = sig32Z << 1;
    uint64_t rem =
        ((sigA - static_cast<uint64_t>(doubleTerm) * static_cast<uint32_t>(sigB >> 32)) << 28) -
        static_cast<uint64_t>(doubleTerm) * (static_cast<uint32_t>(sigB) >> 4);
    uint32_t q = ((static_cast<uint32_t>(rem >> 32) * static_cast<uint64_t>(recip32)) >> 32) + 4;
    uint64_t sigZ = (static_cast<uint64_t>(sig32Z) << 32) + (static_cast<uint64_t>(q) << 4);

    if ((sigZ & 0x1FF) < 4 << 4) {
        q &= ~7;
        sigZ &= ~static_cast<uint64_t>(0x7F);
        doubleTerm = q << 1;
        rem =
            ((rem - static_cast<uint64_t>(doubleTerm) * static_cast<uint32_t>(sigB >> 32)) << 28)
            - static_cast<uint64_t>(doubleTerm) * (static_cast<uint32_t>(sigB) >> 4);

        if (rem & UINT64_C(0x8000000000000000)) {
            sigZ -= 1 << 7;
        } else {
            if (rem) {
                sigZ |= 1;
            }
        }
    }

    return softfloat_roundPackToF64(signZ, expZ, sigZ);
}
