
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

float128_t
f128_sqrt(float128_t const a)
{
    using namespace softfloat::internals::fast_int64;

    uint128 const uA{a};
    bool const signA = is_sign(uA.v64);
    int32_t expA = exp_F128_UI64(uA.v64);
    uint128 sigA{frac_F128_UI64(uA.v64), uA.v0};

    if (0x7FFF == expA) {
        if (0 != (sigA.v64 | sigA.v0)) {
            return float128_t(propagate_NaN(uA, uint128{0, 0}));
        }

        // infinity

        if (!signA) {
            return a;
        }

        // -inf
        softfloat_raiseFlags(softfloat_flag_invalid);
        return float128_t(uint128{defaultNaNF128UI64, defaultNaNF128UI0});
    }

    if (signA) {
        if (0 == (expA | sigA.v64 | sigA.v0)) {
            return a;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return float128_t(uint128{defaultNaNF128UI64, defaultNaNF128UI0});
    }

    if (0 == expA) {
        if (0 == (sigA.v64 | sigA.v0)) {
            // zero
            return a;
        }

        exp32_sig128 const normExpSig = norm_subnormal_F128Sig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    /*
    `sig32Z' is guaranteed to be a lower bound on the square root of
    `sig32A', which makes `sig32Z' also a lower bound on the square root of
    `sigA'.
    */
    int32_t const expZ = ((expA - 0x3FFF) >> 1) + 0x3FFE;
    expA &= 1;
    sigA.v64 |= UINT64_C(0x0001000000000000);
    uint32_t const sig32A = static_cast<uint32_t>(sigA.v64 >> 17);
    uint32_t const recipSqrt32 = approx_recip_sqrt_32_1(static_cast<uint32_t>(expA), sig32A);
    uint32_t sig32Z = (static_cast<uint64_t>(sig32A) * recipSqrt32) >> 32;
    uint128 rem;

    if (expA) {
        sig32Z >>= 1;
        rem = short_shift_left_128(sigA, 12);
    } else {
        rem = short_shift_left_128(sigA, 13);
    }

    uint32_t qs[3];
    qs[2] = sig32Z;
    rem.v64 -= static_cast<uint64_t>(sig32Z) * sig32Z;

    uint32_t const q = (static_cast<uint32_t>(rem.v64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    qs[1] = q;
    uint64_t x64 = static_cast<uint64_t>(sig32Z) << 32;
    uint64_t sig64Z = x64 + (static_cast<uint64_t>(q) << 3);
    x64 += sig64Z;
    uint128 rem_1 = short_shift_left_128(rem, 29);
    uint128 const tmp = mul_64_by_shifted_32_to_128(x64, q);
    uint128 rem_2 = sub(rem_1, tmp);
    uint32_t q_2 = (static_cast<uint32_t>(rem_2.v64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    uint128 const y_0 = short_shift_left_128(rem_2, 29);
    sig64Z <<= 1;

    /* Repeating this loop is a rare occurrence.*/
    for (;;) {
        uint128 const tmp_1 = short_shift_left_128(uint128{0, sig64Z}, 32);
        uint128 const tmp_2 = add(tmp_1, uint128{0, static_cast<uint64_t>(q_2) << 6});
        uint128 const tmp_3 = mul_128_by_32(tmp_2, q_2);
        rem_2 = sub(y_0, tmp_3);

        if (0 == (rem_2.v64 & UINT64_C(0x8000000000000000))) {
            break;
        }

        --q_2;
    }

    qs[0] = q_2;

    uint32_t q_1 = ((static_cast<uint32_t>(rem_2.v64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32) + 2;
    uint64_t sigZExtra = static_cast<uint64_t>(static_cast<uint64_t>(q_1) << 59);
    uint128 const term_1 = short_shift_left_128(uint128{0, qs[1]}, 53);
    uint128 sigZ = add(uint128{static_cast<uint64_t>(qs[2]) << 18, (static_cast<uint64_t>(qs[0]) << 24) + (q_1 >> 5)}, term_1);

    if ((q_1 & 0xF) <= 2) {
        q_1 &= ~3;
        sigZExtra = static_cast<uint64_t>(static_cast<uint64_t>(q_1) << 59);
        uint128 y = short_shift_left_128(sigZ, 6);
        y.v0 |= sigZExtra >> 58;
        uint128 const term_2 = sub(y, uint128{0, q_1});
        uint128 const y_1 = mul_64_by_shifted_32_to_128(term_2.v0, q_1);
        uint128 const term_3 = mul_64_by_shifted_32_to_128(term_2.v64, q_1);
        uint128 const term_4 = add(term_3, uint128{0, y_1.v64});
        uint128 const rem_3 = short_shift_left_128(rem_2, 20);
        uint128 const term_5 = sub(term_4, rem_3);

        /* The concatenation of `term' and `y_1.v0' is now the negative remainder
        (3 words altogether).
        */
        if (0 != (term_5.v64 & UINT64_C(0x8000000000000000))) {
            sigZExtra |= 1;
        } else if (0 != (term_5.v64 | term_5.v0 | y_1.v0)) {
            if (sigZExtra) {
                --sigZExtra;
            } else {
                sigZ = sub(sigZ, uint128{0, 1});
                sigZExtra = ~UINT64_C(0);
            }
        }
    }

    return round_pack_to_F128(false, expZ, sigZ.v64, sigZ.v0, sigZExtra);
}
