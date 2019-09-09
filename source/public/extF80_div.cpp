
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

extFloat80_t
extF80_div(extFloat80_t const a,
           extFloat80_t const b)
{
    using namespace softfloat::internals::fast_int64;

    bool const signA = is_sign(a.signExp);
    int32_t expA = expExtF80UI64(a.signExp);
    uint64_t sigA = a.signif;
    uint16_t const uiB64 = b.signExp;
    uint64_t const uiB0 = b.signif;
    bool const signB = is_sign(uiB64);
    int32_t expB = expExtF80UI64(uiB64);
    uint64_t sigB = uiB0;
    bool const signZ = signA != signB;

    if (0x7FFF == expA) {
        if (sigA & UINT64_C(0x7FFFFFFFFFFFFFFF)) {
            auto const uiZ_1 = propagate_NaN(a.signExp, a.signif, uiB64, uiB0);
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(uiZ_1.v64);
            uZ.signif = uiZ_1.v0;
            return uZ;
        } else if (expB != 0x7FFF) {
            extFloat80_t uZ;
            uZ.signExp = packToExtF80UI64(signZ, 0x7FFF);
            uZ.signif = UINT64_C(0x8000000000000000);
            return uZ;
        } else if (sigB & UINT64_C(0x7FFFFFFFFFFFFFFF)) {
            auto const uiZ_1 = propagate_NaN(a.signExp, a.signif, uiB64, uiB0);
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(uiZ_1.v64);
            uZ.signif = uiZ_1.v0;
            return uZ;
        } else {
            softfloat_raiseFlags(softfloat_flag_invalid);
            extFloat80_t uZ;
            uZ.signExp = defaultNaNExtF80UI64;
            uZ.signif = defaultNaNExtF80UI0;
            return uZ;
        }
    }

    if (0x7FFF == expB) {
        if (sigB & UINT64_C(0x7FFFFFFFFFFFFFFF)) {
            auto const uiZ_1 = propagate_NaN(a.signExp, a.signif, uiB64, uiB0);
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(uiZ_1.v64);
            uZ.signif = uiZ_1.v0;
            return uZ;
        } else {
            extFloat80_t uZ;
            uZ.signExp = packToExtF80UI64(signZ, 0);
            uZ.signif = 0;
            return uZ;
        }
    }

    if (0 == expB) {
        expB = 1;
    }

    if (0 == (sigB & UINT64_C(0x8000000000000000))) {
        if (0 == sigB) {
            if (0 == sigA) {
                softfloat_raiseFlags(softfloat_flag_invalid);
                extFloat80_t uZ;
                uZ.signExp = defaultNaNExtF80UI64;
                uZ.signif = defaultNaNExtF80UI0;
                return uZ;
            }

            softfloat_raiseFlags(softfloat_flag_infinite);
            extFloat80_t uZ;
            uZ.signExp = packToExtF80UI64(signZ, 0x7FFF);
            uZ.signif = UINT64_C(0x8000000000000000);
            return uZ;
        }

        exp32_sig64 const normExpSig = normSubnormalExtF80Sig(sigB);
        expB += normExpSig.exp;
        sigB = normExpSig.sig;
    }

    if (0 == expA) {
        expA = 1;
    }

    if (0 == (UINT64_C(0x8000000000000000) & sigA)) {
        if (0 == sigA) {
            extFloat80_t uZ;
            uZ.signExp = packToExtF80UI64(signZ, 0);
            uZ.signif = 0;
            return uZ;
        }

        exp32_sig64 const normExpSig = normSubnormalExtF80Sig(sigA);
        expA += normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int32_t expZ = expA - expB + 0x3FFF;

    uint128 rem;

    if (sigA < sigB) {
        --expZ;
        rem = shortShiftLeft128(uint128{0, sigA}, 32);
    } else {
        rem = shortShiftLeft128(uint128{0, sigA}, 31);
    }

    uint32_t const recip32 = approxRecip32_1(sigB >> 32);
    uint64_t sigZ = 0;
    uint32_t q;

    for (int ix = 2;;) {
        uint64_t const q64 = static_cast<uint64_t>(static_cast<uint32_t>(rem.v64 >> 2)) * recip32;
        q = (q64 + 0x80000000) >> 32;
        --ix;

        if (ix < 0) {
            break;
        }

        rem = shortShiftLeft128(rem, 29);
        rem = sub(rem, mul64ByShifted32To128(sigB, q));

        if (rem.v64 & UINT64_C(0x8000000000000000)) {
            --q;
            rem = add(rem, uint128{sigB >> 32, sigB << 32});
        }

        sigZ = (sigZ << 29) + q;
    }

    if (((q + 1) & 0x3FFFFF) < 2) {
        rem = shortShiftLeft128(rem, 29);
        rem = sub(rem, mul64ByShifted32To128(sigB, q));
        uint128 const term = shortShiftLeft128(uint128{0, sigB}, 32);

        if (rem.v64 & UINT64_C(0x8000000000000000)) {
            --q;
            rem = add(rem, term);
        } else if (le(term, rem)) {
            ++q;
            rem = sub(rem, term);
        }

        if (rem.v64 | rem.v0) {
            q |= 1;
        }
    }

    return roundPackToExtF80(signZ, expZ, (sigZ << 6) + (q >> 23), static_cast<uint64_t>(static_cast<uint64_t>(q) << 41), extF80_roundingPrecision);
}
