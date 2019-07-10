
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

#include "target.hpp"

extFloat80_t
extF80_mul(extFloat80_t const a,
           extFloat80_t const b)
{
    using namespace softfloat::internals;

    bool const signA = signExtF80UI64(a.signExp);
    int32_t const expA = expExtF80UI64(a.signExp);
    uint64_t sigA = a.signif;
    bool const signB = signExtF80UI64(b.signExp);
    int32_t const expB = expExtF80UI64(b.signExp);
    uint64_t sigB = b.signif;
    bool const signZ = signA != signB;

    if (0x7FFF == expA) {
        extFloat80_t uZ_1;

        if (0 != (UINT64_C(0x7FFFFFFFFFFFFFFF) & sigA) || (expB == 0x7FFF && 0 != (UINT64_C(0x7FFFFFFFFFFFFFFF) & sigB))) {
            auto const uiZ_1 = softfloat_propagateNaNExtF80UI(a.signExp, a.signif, b.signExp, b.signif);
            uZ_1.signExp = static_cast<uint16_t>(uiZ_1.v64);
            uZ_1.signif = uiZ_1.v0;
        } else if (0 == (expB | sigB)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            uZ_1.signExp = defaultNaNExtF80UI64;
            uZ_1.signif = defaultNaNExtF80UI0;
        } else {
            uZ_1.signExp = packToExtF80UI64(signZ, 0x7FFF);
            uZ_1.signif = UINT64_C(0x8000000000000000);
        }

        return uZ_1;
    } else if (0x7FFF == expB) {
        extFloat80_t uZ_1;

        if (0 != (UINT64_C(0x7FFFFFFFFFFFFFFF) & sigB)) {
            auto const uiZ_1 = softfloat_propagateNaNExtF80UI(a.signExp, a.signif, b.signExp, b.signif);
            uZ_1.signExp = static_cast<uint16_t>(uiZ_1.v64);
            uZ_1.signif = uiZ_1.v0;
        } else if (0 == (expA | sigA)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            uZ_1.signExp = defaultNaNExtF80UI64;
            uZ_1.signif = defaultNaNExtF80UI0;
        } else {
            uZ_1.signExp = packToExtF80UI64(signZ, 0x7FFF);
            uZ_1.signif = UINT64_C(0x8000000000000000);
        }

        return uZ_1;
    } else {
        auto expA_1 =  0 == expA ? 1 : expA;

        if (0 == (sigA & UINT64_C(0x8000000000000000))) {
            if (0 == sigA) {
                extFloat80_t uZ_1;
                uZ_1.signExp = packToExtF80UI64(signZ, 0);
                uZ_1.signif = 0;
                return uZ_1;
            } else {
                exp32_sig64 const normExpSig = softfloat_normSubnormalExtF80Sig(sigA);
                expA_1 += normExpSig.exp;
                sigA = normExpSig.sig;
            }
        }

        auto expB_1 = 0 == expB ? 1 : expB;

        if (0 == (UINT64_C(0x8000000000000000) & sigB)) {
            if (0 == sigB) {
                extFloat80_t uZ_1;
                uZ_1.signExp = packToExtF80UI64(signZ, 0);
                uZ_1.signif = 0;
                return uZ_1;
            } else {
                exp32_sig64 const normExpSig = softfloat_normSubnormalExtF80Sig(sigB);
                expB_1 += normExpSig.exp;
                sigB = normExpSig.sig;
            }
        }

        int32_t const expZ = expA_1 + expB_1 - 0x3FFE;
        uint128 const sig128Z = softfloat_mul64To128(sigA, sigB);

        if (sig128Z.v64 < UINT64_C(0x8000000000000000)) {
            auto const sig128Z_1 = softfloat_add128(sig128Z.v64, sig128Z.v0, sig128Z.v64, sig128Z.v0);
            return
                softfloat_roundPackToExtF80(signZ, expZ - 1, sig128Z_1.v64, sig128Z_1.v0, extF80_roundingPrecision);
        } else {
            return
                softfloat_roundPackToExtF80(signZ, expZ, sig128Z.v64, sig128Z.v0, extF80_roundingPrecision);
        }
    }
}

