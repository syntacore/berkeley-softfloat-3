
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

extFloat80_t
extF80_rem(extFloat80_t const a,
           extFloat80_t const b)
{
    using namespace softfloat::internals::fast_int64;

    bool const signA = is_sign(a.signExp);
    int32_t expA = expExtF80UI64(a.signExp);
    uint64_t sigA = a.signif;
    int32_t expB = expExtF80UI64(b.signExp);
    uint64_t sigB = b.signif;

    if (0x7FFF == expA) {
        if (0 != (sigA & UINT64_C(0x7FFFFFFFFFFFFFFF)) || (0x7FFF == expB && 0 != (sigB & UINT64_C(0x7FFFFFFFFFFFFFFF)))) {
            uint128 const uiZ = propagate_NaN(a.signExp, a.signif, b.signExp, b.signif);
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(uiZ.v64);
            uZ.signif = uiZ.v0;
            return uZ;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        extFloat80_t uZ;
        uZ.signExp = defaultNaNExtF80UI64;
        uZ.signif = defaultNaNExtF80UI0;
        return uZ;
    }

    if (0x7FFF == expB) {
        if (0 != (sigB & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
            uint128 const uiZ = propagate_NaN(a.signExp, a.signif, b.signExp, b.signif);
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(uiZ.v64);
            uZ.signif = uiZ.v0;
            return uZ;
        }

        /*
        Argument b is an infinity.  Doubling `expB' is an easy way to ensure
        that `expDiff' later is less than -1, which will result in returning
        a canonicalized version of argument a.
        */
        expB += expB;
    }

    if (0 == expB) {
        expB = 1;
    }

    if (0 == (sigB & UINT64_C(0x8000000000000000))) {
        if (0 == sigB) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            extFloat80_t uZ;
            uZ.signExp = defaultNaNExtF80UI64;
            uZ.signif = defaultNaNExtF80UI0;
            return uZ;
        }

        exp32_sig64 const normExpSig = softfloat_normSubnormalExtF80Sig(sigB);
        expB += normExpSig.exp;
        sigB = normExpSig.sig;
    }

    if (0 == expA) {
        expA = 1;
    }

    if (0 == (sigA & UINT64_C(0x8000000000000000))) {
        if (0 == sigA) {
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(packToExtF80UI64(signA, 0));
            uZ.signif = sigA >> 1;
            return uZ;
        }

        exp32_sig64 const normExpSig = softfloat_normSubnormalExtF80Sig(sigA);
        expA += normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int32_t expDiff = expA - expB;

    if (expDiff < -1) {
        if (expA < 1) {
            sigA >>= 1 - expA;
            expA = 0;
        }

        extFloat80_t uZ;
        uZ.signExp = static_cast<uint16_t>(packToExtF80UI64(signA, static_cast<uint16_t>(expA)));
        uZ.signif = sigA;
        return uZ;
    }

    uint128 rem = softfloat_shortShiftLeft128(uint128{0, sigA}, 32);
    uint128 shiftedSigB = softfloat_shortShiftLeft128(uint128{0, sigB}, 32);

    uint128 altRem;
    uint32_t q;

    if (expDiff < 1) {
        if (0 != expDiff) {
            --expB;
            shiftedSigB = softfloat_shortShiftLeft128(uint128{0, sigB}, 33);
            q = 0;
        } else {
            q = 0u + !!(sigB <= sigA);

            if (0 != q) {
                rem = softfloat_sub128(rem, shiftedSigB);
            }
        }
    } else {
        uint32_t const recip32 = softfloat_approxRecip32_1(sigB >> 32);
        expDiff -= 30;

        uint64_t q64;

        for (;;) {
            q64 = static_cast<uint64_t>(static_cast<uint32_t>(rem.v64 >> 2)) * recip32;

            if (expDiff < 0) {
                break;
            }

            q = (q64 + 0x80000000) >> 32;
            rem = softfloat_shortShiftLeft128(rem, 29);
            uint128 term = softfloat_mul64ByShifted32To128(sigB, q);
            rem = softfloat_sub128(rem, term);

            if (0 != (rem.v64 & UINT64_C(0x8000000000000000))) {
                rem = softfloat_add128(rem, shiftedSigB);
            }

            expDiff -= 29;
        }

        /* `expDiff' cannot be less than -29 here.*/
        assert(-29 <= expDiff);
        q = static_cast<uint32_t>(q64 >> 32) >> (~expDiff & 31);
        /**
        @todo Warning   C4244   '=': conversion from 'int' to 'uint8_t', possible loss of data
        */
        rem = softfloat_shortShiftLeft128(rem, uint8_t(expDiff + 30));
        rem = softfloat_sub128(rem, softfloat_mul64ByShifted32To128(sigB, q));

        if (0 != (rem.v64 & UINT64_C(0x8000000000000000))) {
            altRem = softfloat_add128(rem, shiftedSigB);
            uint128 const meanRem = softfloat_add128(rem, altRem);

            if (0 != (meanRem.v64 & UINT64_C(0x8000000000000000)) || (!(meanRem.v64 | meanRem.v0) && (q & 1))) {
                rem = altRem;
            }

            bool signRem = signA;

            if (0 != (rem.v64 & UINT64_C(0x8000000000000000))) {
                signRem = !signRem;
                rem = softfloat_sub128(uint128{0, 0}, rem);
            }

            return softfloat_normRoundPackToExtF80(signRem, expB + 32, rem.v64, rem.v0, 80);
        }
    }

    do {
        altRem = rem;
        ++q;
        rem = softfloat_sub128(rem, shiftedSigB);
    } while (0 == (rem.v64 & UINT64_C(0x8000000000000000)));

    uint128 const meanRem = softfloat_add128(rem, altRem);

    if (
        0 != (meanRem.v64 & UINT64_C(0x8000000000000000)) ||
        (0 == (meanRem.v64 | meanRem.v0) && 0 != (q & 1))
    ) {
        rem = altRem;
    }

    bool signRem = signA;

    if (0 != (rem.v64 & UINT64_C(0x8000000000000000))) {
        signRem = !signRem;
        rem = softfloat_sub128(uint128{0, 0}, rem);
    }

    return softfloat_normRoundPackToExtF80(signRem, expB + 32, rem.v64, rem.v0, 80);
}
