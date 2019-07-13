
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

extFloat80_t
extF80_sqrt(extFloat80_t const a)
{
    using namespace softfloat::internals;

    bool const signA = is_sign(a.signExp);
    int32_t expA = expExtF80UI64(a.signExp);
    uint64_t sigA = a.signif;

    if (0x7FFF == expA) {
        if (0 != (UINT64_C(0x7FFFFFFFFFFFFFFF) & sigA)) {
            uint128 const uiZ = softfloat_propagateNaNExtF80UI(a.signExp, a.signif, 0, 0);
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(uiZ.v64);
            uZ.signif = uiZ.v0;
            return uZ;
        }

        if (!signA) {
            return a;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        extFloat80_t uZ;
        uZ.signExp = defaultNaNExtF80UI64;
        uZ.signif = defaultNaNExtF80UI0;
        return uZ;
    }

    if (signA) {
        if (0 == sigA) {
            extFloat80_t uZ;
            uZ.signExp = packToExtF80UI64(signA, 0);
            uZ.signif = 0;
            return uZ;
        } else {
            softfloat_raiseFlags(softfloat_flag_invalid);
            extFloat80_t uZ;
            uZ.signExp = defaultNaNExtF80UI64;
            uZ.signif = defaultNaNExtF80UI0;
            return uZ;
        }
    }

    if (0 == expA) {
        expA = 1;
    }

    if (0 == (UINT64_C(0x8000000000000000) & sigA)) {
        if (!sigA) {
            extFloat80_t uZ;
            uZ.signExp = packToExtF80UI64(signA, 0);
            uZ.signif = 0;
            return uZ;
        }

        exp32_sig64 const normExpSig = softfloat_normSubnormalExtF80Sig(sigA);
        expA += normExpSig.exp;
        sigA = normExpSig.sig;
    }

    /*
    `sig32Z' is guaranteed to be a lower bound on the square root of
    `sig32A', which makes `sig32Z' also a lower bound on the square root of
    `sigA'.)
    */
    int32_t const expZ = ((expA - 0x3FFF) >> 1) + 0x3FFF;
    expA &= 1;
    uint32_t const sig32A = sigA >> 32;
    uint32_t const recipSqrt32 = softfloat_approxRecipSqrt32_1(static_cast<uint32_t>(expA), sig32A);
    uint32_t sig32Z = (static_cast<uint64_t>(sig32A) * recipSqrt32) >> 32;

    uint128 rem;

    if (0 != expA) {
        sig32Z >>= 1;
        rem = softfloat_shortShiftLeft128(0, sigA, 61);
    } else {
        rem = softfloat_shortShiftLeft128(0, sigA, 62);
    }

    rem.v64 -= static_cast<uint64_t>(sig32Z) * sig32Z;

    uint64_t q = (static_cast<uint32_t>(rem.v64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    uint64_t sigZ = (static_cast<uint64_t>(sig32Z) << 32) + (q << 3);
    uint64_t x64 = (static_cast<uint64_t>(sig32Z) << 32) + sigZ;
    /**
    @todo Warning   C4242   'function': conversion from 'int64_t' to 'int32_t', possible loss of data
    */
    uint128 term = softfloat_mul64ByShifted32To128(x64, static_cast<uint32_t>(q));
    rem = softfloat_shortShiftLeft128(rem, 29);
    rem = softfloat_sub128(rem, term);

    q = ((static_cast<uint32_t>(rem.v64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32) + 2;
    x64 = sigZ;
    sigZ = (sigZ << 1) + (q >> 25);
    uint64_t sigZExtra = static_cast<uint64_t>(q << 39);

    if ((q & 0xFFFFFF) <= 2) {
        q &= ~static_cast<uint64_t>(0xFFFF);
        sigZExtra = static_cast<uint64_t>(q << 39);
        /**
        @todo Warning   C4242   'function': conversion from 'int64_t' to 'int32_t', possible loss of data
        */
        term = softfloat_mul64ByShifted32To128(x64 + (q >> 27), static_cast<uint32_t>(q));
        x64 = static_cast<uint32_t>(q << 5) * static_cast<uint64_t>(static_cast<uint32_t>(q));
        term = softfloat_add128(term.v64, term.v0, 0, x64);
        rem = softfloat_shortShiftLeft128(rem, 28);
        rem = softfloat_sub128(rem, term);

        if (0 != (UINT64_C(0x8000000000000000) & rem.v64)) {
            if (0 == sigZExtra) {
                --sigZ;
            }

            --sigZExtra;
        } else if (0 != (rem.v64 | rem.v0)) {
            sigZExtra |= 1;
        }
    }

    return softfloat_roundPackToExtF80(0, expZ, sigZ, sigZExtra, extF80_roundingPrecision);
}

