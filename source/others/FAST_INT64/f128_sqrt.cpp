
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

#ifndef SOFTFLOAT_FAST_INT64
#error Fast int64_t operations only
#endif

float128_t
f128_sqrt(float128_t const a)
{
    using namespace softfloat::internals;
    uint128 const uA{a};
    bool const signA = is_sign(uA.v64);
    int32_t expA = expF128UI64(uA.v64);
    uint128 sigA{fracF128UI64(uA.v64), uA.v0};

    if (0x7FFF == expA) {
        if (0 != (sigA.v64 | sigA.v0)) {
            return float128_t(softfloat_propagateNaNF128UI(uA.v64, uA.v0, 0, 0));
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

        exp32_sig128 const normExpSig = softfloat_normSubnormalF128Sig(sigA.v64, sigA.v0);
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
    uint32_t const recipSqrt32 = softfloat_approxRecipSqrt32_1(static_cast<uint32_t>(expA), sig32A);
    uint32_t sig32Z = (static_cast<uint64_t>(sig32A) * recipSqrt32) >> 32;
    uint128 rem;

    if (expA) {
        sig32Z >>= 1;
        rem = softfloat_shortShiftLeft128(sigA, 12);
    } else {
        rem = softfloat_shortShiftLeft128(sigA, 13);
    }

    uint32_t qs[3];
    qs[2] = sig32Z;
    rem.v64 -= static_cast<uint64_t>(sig32Z) * sig32Z;

    uint32_t const q = (static_cast<uint32_t>(rem.v64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    qs[1] = q;
    uint64_t x64 = static_cast<uint64_t>(sig32Z) << 32;
    uint64_t sig64Z = x64 + (static_cast<uint64_t>(q) << 3);
    x64 += sig64Z;
    uint128 rem_1 = softfloat_shortShiftLeft128(rem, 29);
    uint128 const tmp = softfloat_mul64ByShifted32To128(x64, q);
    uint128 rem_2 = softfloat_sub128(rem_1.v64, rem_1.v0, tmp.v64, tmp.v0);
    uint32_t q_2 = (static_cast<uint32_t>(rem_2.v64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    uint128 const y_0 = softfloat_shortShiftLeft128(rem_2, 29);
    sig64Z <<= 1;

    /* Repeating this loop is a rare occurrence.*/
    for (;;) {
        uint128 const tmp_1 = softfloat_shortShiftLeft128(0, sig64Z, 32);
        uint128 const tmp_2 = softfloat_add128(tmp_1.v64, tmp_1.v0, 0, static_cast<uint64_t>(q_2) << 6);
        uint128 const tmp_3 = softfloat_mul128By32(tmp_2.v64, tmp_2.v0, q_2);
        rem_2 = softfloat_sub128(y_0.v64, y_0.v0, tmp_3.v64, tmp_3.v0);

        if (0 == (rem_2.v64 & UINT64_C(0x8000000000000000))) {
            break;
        }

        --q_2;
    }

    qs[0] = q_2;

    uint32_t q_1 = ((static_cast<uint32_t>(rem_2.v64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32) + 2;
    uint64_t sigZExtra = static_cast<uint64_t>(static_cast<uint64_t>(q_1) << 59);
    uint128 const term_1 = softfloat_shortShiftLeft128(0, qs[1], 53);
    uint128 sigZ = softfloat_add128(static_cast<uint64_t>(qs[2]) << 18, (static_cast<uint64_t>(qs[0]) << 24) + (q_1 >> 5), term_1.v64, term_1.v0);

    if ((q_1 & 0xF) <= 2) {
        q_1 &= ~3;
        sigZExtra = static_cast<uint64_t>(static_cast<uint64_t>(q_1) << 59);
        uint128 y = softfloat_shortShiftLeft128(sigZ, 6);
        y.v0 |= sigZExtra >> 58;
        uint128 const term_2 = softfloat_sub128(y.v64, y.v0, 0, q_1);
        uint128 const y_1 = softfloat_mul64ByShifted32To128(term_2.v0, q_1);
        uint128 const term_3 = softfloat_mul64ByShifted32To128(term_2.v64, q_1);
        uint128 const term_4 = softfloat_add128(term_3.v64, term_3.v0, 0, y_1.v64);
        uint128 const rem_3 = softfloat_shortShiftLeft128(rem_2, 20);
        uint128 const term_5 = softfloat_sub128(term_4.v64, term_4.v0, rem_3.v64, rem_3.v0);

        /* The concatenation of `term' and `y_1.v0' is now the negative remainder
        (3 words altogether).
        */
        if (0 != (term_5.v64 & UINT64_C(0x8000000000000000))) {
            sigZExtra |= 1;
        } else if (0 != (term_5.v64 | term_5.v0 | y_1.v0)) {
            if (sigZExtra) {
                --sigZExtra;
            } else {
                sigZ = softfloat_sub128(sigZ.v64, sigZ.v0, 0, 1);
                sigZExtra = ~UINT64_C(0);
            }
        }
    }

    return softfloat_roundPackToF128(false, expZ, sigZ.v64, sigZ.v0, sigZExtra);
}
