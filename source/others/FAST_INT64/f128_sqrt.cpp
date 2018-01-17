
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

#include "softfloat/functions.h"

#include "internals.hpp"
#include "specialize.hpp"

float128_t
f128_sqrt(float128_t a)
{
    using namespace softfloat;
    ui128_f128 uA;
    uA.f = a;
    uint64_t const uiA64 = uA.ui.v64;
    uint64_t const uiA0 = uA.ui.v0;
    bool const signA = signF128UI64(uiA64);
    int32_t expA = expF128UI64(uiA64);
    uint128 sigA;
    sigA.v64 = fracF128UI64(uiA64);
    sigA.v0 = uiA0;

    if (expA == 0x7FFF) {
        if (sigA.v64 | sigA.v0) {
            return u_as_f_128(softfloat_propagateNaNF128UI(uiA64, uiA0, 0, 0));
        }

        if (!signA) {
            return a;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        uint128 uiZ;
        uiZ.v64 = defaultNaNF128UI64;
        uiZ.v0 = defaultNaNF128UI0;
        return u_as_f_128(uiZ);
    }

    if (signA) {
        if (!(expA | sigA.v64 | sigA.v0)) {
            return a;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        uint128 uiZ;
        uiZ.v64 = defaultNaNF128UI64;
        uiZ.v0 = defaultNaNF128UI0;
        return u_as_f_128(uiZ);
    } else {
        if (!expA) {
            if (!(sigA.v64 | sigA.v0)) {
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
        uint32_t sig32Z = ((uint64_t)sig32A * recipSqrt32) >> 32;
        uint128 rem;

        if (expA) {
            sig32Z >>= 1;
            rem = softfloat_shortShiftLeft128(sigA.v64, sigA.v0, 12);
        } else {
            rem = softfloat_shortShiftLeft128(sigA.v64, sigA.v0, 13);
        }

        uint32_t qs[3];
        qs[2] = sig32Z;
        rem.v64 -= (uint64_t)sig32Z * sig32Z;

        uint32_t q = ((uint32_t)(rem.v64 >> 2) * (uint64_t)recipSqrt32) >> 32;
        qs[1] = q;
        uint64_t x64 = (uint64_t)sig32Z << 32;
        uint64_t sig64Z = x64 + ((uint64_t)q << 3);
        x64 += sig64Z;
        rem = softfloat_shortShiftLeft128(rem.v64, rem.v0, 29);
        uint128 term = softfloat_mul64ByShifted32To128(x64, q);
        rem = softfloat_sub128(rem.v64, rem.v0, term.v64, term.v0);

        q = ((uint32_t)(rem.v64 >> 2) * (uint64_t)recipSqrt32) >> 32;
        uint128 y = softfloat_shortShiftLeft128(rem.v64, rem.v0, 29);
        sig64Z <<= 1;

        /* Repeating this loop is a rare occurrence.*/
        for (;;) {
            term = softfloat_shortShiftLeft128(0, sig64Z, 32);
            term = softfloat_add128(term.v64, term.v0, 0, (uint64_t)q << 6);
            term = softfloat_mul128By32(term.v64, term.v0, q);
            rem = softfloat_sub128(y.v64, y.v0, term.v64, term.v0);

            if (!(rem.v64 & UINT64_C(0x8000000000000000))) {
                break;
            }

            --q;
        }

        qs[0] = q;

        q = (((uint32_t)(rem.v64 >> 2) * (uint64_t)recipSqrt32) >> 32) + 2;
        uint64_t sigZExtra = (uint64_t)((uint64_t)q << 59);
        term = softfloat_shortShiftLeft128(0, qs[1], 53);
        uint128 sigZ = softfloat_add128((uint64_t)qs[2] << 18, ((uint64_t)qs[0] << 24) + (q >> 5), term.v64, term.v0);

        if ((q & 0xF) <= 2) {
            q &= ~3;
            sigZExtra = (uint64_t)((uint64_t)q << 59);
            y = softfloat_shortShiftLeft128(sigZ.v64, sigZ.v0, 6);
            y.v0 |= sigZExtra >> 58;
            term = softfloat_sub128(y.v64, y.v0, 0, q);
            y = softfloat_mul64ByShifted32To128(term.v0, q);
            term = softfloat_mul64ByShifted32To128(term.v64, q);
            term = softfloat_add128(term.v64, term.v0, 0, y.v64);
            rem = softfloat_shortShiftLeft128(rem.v64, rem.v0, 20);
            term = softfloat_sub128(term.v64, term.v0, rem.v64, rem.v0);

            /* The concatenation of `term' and `y.v0' is now the negative remainder
            (3 words altogether).
            */
            if (term.v64 & UINT64_C(0x8000000000000000)) {
                sigZExtra |= 1;
            } else {
                if (term.v64 | term.v0 | y.v0) {
                    if (sigZExtra) {
                        --sigZExtra;
                    } else {
                        sigZ = softfloat_sub128(sigZ.v64, sigZ.v0, 0, 1);
                        sigZExtra = ~UINT64_C(0);
                    }
                }
            }
        }

        return softfloat_roundPackToF128(0, expZ, sigZ.v64, sigZ.v0, sigZExtra);
    }
}

