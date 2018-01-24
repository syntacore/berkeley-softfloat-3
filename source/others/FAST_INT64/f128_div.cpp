
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

#include "softfloat/functions.h"

#include "internals.hpp"
#include "target.hpp"

float128_t
f128_div(float128_t a,
         float128_t b)
{
    using namespace softfloat::internals;
    uint64_t const uiA64 = f_as_u_128(a).v64;
    uint64_t const uiA0 = f_as_u_128(a).v0;
    bool const signA = signF128UI64(uiA64);
    int32_t expA = expF128UI64(uiA64);

    uint128 sigA;
    sigA.v64 = fracF128UI64(uiA64);
    sigA.v0 = uiA0;

    uint64_t const uiB64 = f_as_u_128(b).v64;
    uint64_t const uiB0 = f_as_u_128(b).v0;
    bool const signB = signF128UI64(uiB64);
    int32_t expB = expF128UI64(uiB64);

    uint128 sigB;
    sigB.v64 = fracF128UI64(uiB64);
    sigB.v0 = uiB0;

    bool const signZ = signA ^ signB;

    if (expA == 0x7FFF) {
        if (sigA.v64 | sigA.v0) {
            uint128 const uiZ = softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0);
            return u_as_f_128(uiZ);
        }

        if (expB == 0x7FFF) {
            if (sigB.v64 | sigB.v0) {
                uint128 const uiZ = softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0);
                return u_as_f_128(uiZ);
            }

            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_128(uint128{defaultNaNF128UI64, defaultNaNF128UI0});
        }

        return u_as_f_128(uint128{packToF128UI64(signZ, 0x7FFF, 0), 0});
    }

    if (expB == 0x7FFF) {
        if (sigB.v64 | sigB.v0) {
            return u_as_f_128(softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0));
        }

        return u_as_f_128(uint128{packToF128UI64(signZ, 0, 0), 0});
    }

    if (!expB) {
        if (!(sigB.v64 | sigB.v0)) {
            if (!(expA | sigA.v64 | sigA.v0)) {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return u_as_f_128(uint128{defaultNaNF128UI64, defaultNaNF128UI0});
            }

            softfloat_raiseFlags(softfloat_flag_infinite);
            return u_as_f_128(uint128{packToF128UI64(signZ, 0x7FFF, 0), 0});
        }

        exp32_sig128 const normExpSig = softfloat_normSubnormalF128Sig(sigB.v64, sigB.v0);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    if (!expA) {
        if (!(sigA.v64 | sigA.v0)) {
            return u_as_f_128(uint128{packToF128UI64(signZ, 0, 0), 0});
        }

        exp32_sig128 const normExpSig = softfloat_normSubnormalF128Sig(sigA.v64, sigA.v0);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int32_t expZ = expA - expB + 0x3FFE;
    sigA.v64 |= UINT64_C(0x0001000000000000);
    sigB.v64 |= UINT64_C(0x0001000000000000);
    uint128 rem = sigA;

    if (softfloat_lt128(sigA.v64, sigA.v0, sigB.v64, sigB.v0)) {
        --expZ;
        rem = softfloat_add128(sigA.v64, sigA.v0, sigA.v64, sigA.v0);
    }

    uint32_t const recip32 = softfloat_approxRecip32_1(static_cast<uint32_t>(sigB.v64 >> 17));
    int ix = 3;

    uint32_t q;
    uint32_t qs[3];

    for (;;) {
        auto const q64 = static_cast<uint64_t>(static_cast<uint32_t>(rem.v64 >> 19)) * recip32;
        q = (q64 + 0x80000000) >> 32;
        --ix;

        if (ix < 0) {
            break;
        }

        rem = softfloat_shortShiftLeft128(rem.v64, rem.v0, 29);
        uint128 const term = softfloat_mul128By32(sigB.v64, sigB.v0, q);
        rem = softfloat_sub128(rem.v64, rem.v0, term.v64, term.v0);

        if (rem.v64 & UINT64_C(0x8000000000000000)) {
            --q;
            rem = softfloat_add128(rem.v64, rem.v0, sigB.v64, sigB.v0);
        }

        qs[ix] = q;
    }

    if (((q + 1) & 7) < 2) {
        rem = softfloat_shortShiftLeft128(rem.v64, rem.v0, 29);
        uint128 const term = softfloat_mul128By32(sigB.v64, sigB.v0, q);
        rem = softfloat_sub128(rem.v64, rem.v0, term.v64, term.v0);

        if (rem.v64 & UINT64_C(0x8000000000000000)) {
            --q;
            rem = softfloat_add128(rem.v64, rem.v0, sigB.v64, sigB.v0);
        } else if (softfloat_le128(sigB.v64, sigB.v0, rem.v64, rem.v0)) {
            ++q;
            rem = softfloat_sub128(rem.v64, rem.v0, sigB.v64, sigB.v0);
        }

        if (rem.v64 | rem.v0) {
            q |= 1;
        }
    }

    uint64_t const sigZExtra = static_cast<uint64_t>(q) << 60;
    uint128 const term = softfloat_shortShiftLeft128(0, qs[1], 54);
    uint128 const sigZ =
        softfloat_add128(
            static_cast<uint64_t>(qs[2]) << 19, (static_cast<uint64_t>(qs[0]) << 25) + (q >> 4),
            term.v64,
            term.v0);
    return softfloat_roundPackToF128(signZ, expZ, sigZ.v64, sigZ.v0, sigZExtra);
}
