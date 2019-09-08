
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

#ifdef _MSC_VER
#pragma warning( disable : 5030)
#endif

extFloat80_t
extF80_roundToInt(extFloat80_t const a,
                  softfloat_round_mode const roundingMode,
                  bool const exact)
{
    using namespace softfloat::internals::fast_int64;

    uint16_t const uiA64 = a.signExp;
    uint16_t const signUI64 = static_cast<uint16_t>(uiA64 & packToExtF80UI64(1, 0));
    int32_t exp = expExtF80UI64(uiA64);
    uint64_t sigA = a.signif;

    if (0 == (sigA & UINT64_C(0x8000000000000000)) && 0x7FFF != exp) {
        if (0 == sigA) {
            extFloat80_t uZ;
            uZ.signExp = signUI64;
            uZ.signif = 0;
            return uZ;
        }

        exp32_sig64 const normExpSig = softfloat_normSubnormalExtF80Sig(sigA);
        exp += normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (0x403E <= exp) {
        uint64_t sigZ;

        if (exp == 0x7FFF) {
            if (sigA & UINT64_C(0x7FFFFFFFFFFFFFFF)) {
                uint128 const uiZ = propagate_NaN(uiA64, sigA, 0, 0);
                extFloat80_t uZ;
                uZ.signExp = static_cast<uint16_t>(uiZ.v64);
                uZ.signif = uiZ.v0;
                return uZ;
            }

            sigZ = UINT64_C(0x8000000000000000);
        } else {
            sigZ = sigA;
        }

        extFloat80_t uZ;
        uZ.signExp = static_cast<uint16_t>(signUI64 | exp);
        uZ.signif = sigZ;
        return uZ;
    }

    if (exp <= 0x3FFE) {
        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        /**
        @bug  warning: enumeration value ‘softfloat_round_minMag’ not handled in switch
        */
        switch (roundingMode) {
        case softfloat_round_near_even:
            if (0 == (sigA & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
                break;
            }

            [[fallthrough]];

        case softfloat_round_near_maxMag:
            if (0x3FFE == exp) {
                extFloat80_t uZ;
                uZ.signExp = signUI64 | 0x3FFFu;
                uZ.signif = UINT64_C(0x8000000000000000);
                return uZ;
            }

            break;

        case softfloat_round_min:
            if (0 != signUI64) {
                extFloat80_t uZ;
                uZ.signExp = signUI64 | 0x3FFFu;
                uZ.signif = UINT64_C(0x8000000000000000);
                return uZ;
            }

            break;

        case softfloat_round_max:
            if (0 == signUI64) {
                extFloat80_t uZ;
                uZ.signExp = signUI64 | 0x3FFFu;
                uZ.signif = UINT64_C(0x8000000000000000);
                return uZ;
            }

            break;
        }

        extFloat80_t uZ;
        uZ.signExp = signUI64;
        uZ.signif = 0;
        return uZ;
    }

    uint16_t uiZ64 = static_cast<uint16_t>(signUI64 | exp);
    uint64_t const lastBitMask = static_cast<uint64_t>(1) << (0x403E - exp);
    uint64_t const roundBitsMask = lastBitMask - 1;
    uint64_t sigZ = sigA;

    if (softfloat_round_near_maxMag == roundingMode) {
        sigZ += lastBitMask >> 1;
    } else if (softfloat_round_near_even == roundingMode) {
        sigZ += lastBitMask >> 1;

        if (!(sigZ & roundBitsMask)) {
            sigZ &= ~lastBitMask;
        }
    } else if (softfloat_round_minMag != roundingMode) {
        if ((0 != signUI64) != (softfloat_round_max == roundingMode)) {
            sigZ += roundBitsMask;
        }
    }

    sigZ &= ~roundBitsMask;

    if (!sigZ) {
        ++uiZ64;
        sigZ = UINT64_C(0x8000000000000000);
    }

    if (exact && sigZ != sigA) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    extFloat80_t uZ;
    uZ.signExp = uiZ64;
    uZ.signif = sigZ;
    return uZ;
}

