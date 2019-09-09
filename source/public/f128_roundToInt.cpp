
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

float128_t
f128_roundToInt(float128_t const a,
                uint8_t const roundingMode,
                bool const exact)
{
    using namespace softfloat::internals::fast_int64;

    uint128 uiZ;

    uint128 const uA{a};
    int32_t exp = exp_F128_UI64(uA.v64);

    if (0x402F <= exp) {
        if (0x406F <= exp) {
            if (0x7FFF == exp && 0 != (frac_F128_UI64(uA.v64) | uA.v0)) {
                return float128_t(propagate_NaN(uA, uint128{0, 0}));
            }

            return a;
        }

        uint64_t const lastBitMask = static_cast<uint64_t>(2) << (0x406E - exp);
        uint64_t const roundBitsMask = lastBitMask - 1;
        uiZ = uA;
        bool const roundNearEven = (roundingMode == softfloat_round_near_even);

        if (roundNearEven || softfloat_round_near_maxMag == roundingMode) {
            if (0x402F == exp) {
                if (UINT64_C(0x8000000000000000) <= uiZ.v0) {
                    ++uiZ.v64;

                    if (roundNearEven && UINT64_C(0x8000000000000000) == uiZ.v0) {
                        uiZ.v64 &= ~1;
                    }
                }
            } else {
                uiZ = add(uiZ, uint128{0, lastBitMask >> 1});

                if (roundNearEven && 0 == (uiZ.v0 & roundBitsMask)) {
                    uiZ.v0 &= ~lastBitMask;
                }
            }
        } else if (softfloat_round_minMag != roundingMode) {
            if (is_sign(uiZ.v64) != (softfloat_round_max == roundingMode)) {
                uiZ = add(uiZ, uint128{0, roundBitsMask});
            }
        }

        uiZ.v0 &= ~roundBitsMask;
    } else {
        if (exp < 0x3FFF) {
            if (0 == ((uA.v64 & UINT64_C(0x7FFFFFFFFFFFFFFF)) | uA.v0)) {
                return a;
            }

            if (exact) {
                softfloat_raiseFlags(softfloat_flag_inexact);
            }

            uiZ.v64 = uA.v64 & pack_to_F128_UI64(1, 0, 0);
            uiZ.v0 = 0;

            /**
            @bug  warning: enumeration value ‘softfloat_round_minMag’ not handled in switch
            */
            switch (roundingMode) {
            case softfloat_round_near_even:
                if (0 == (frac_F128_UI64(uA.v64) | uA.v0)) {
                    break;
                }

                [[fallthrough]];

            case softfloat_round_near_maxMag:
                if (0x3FFE == exp) {
                    uiZ.v64 |= pack_to_F128_UI64(0, 0x3FFF, 0);
                }

                break;

            case softfloat_round_min:
                if (0 != uiZ.v64) {
                    uiZ.v64 = pack_to_F128_UI64(1, 0x3FFF, 0);
                }

                break;

            case softfloat_round_max:
                if (0 == uiZ.v64) {
                    uiZ.v64 = pack_to_F128_UI64(0, 0x3FFF, 0);
                }

                break;
            }

            return static_cast<float128_t>(uiZ);
        }

        uiZ.v64 = uA.v64;
        uiZ.v0 = 0;
        uint64_t const lastBitMask = static_cast<uint64_t>(1) << (0x402F - exp);
        uint64_t const roundBitsMask = lastBitMask - 1;

        if (roundingMode == softfloat_round_near_maxMag) {
            uiZ.v64 += lastBitMask >> 1;
        } else if (roundingMode == softfloat_round_near_even) {
            uiZ.v64 += lastBitMask >> 1;

            if (0 == ((uiZ.v64 & roundBitsMask) | uA.v0)) {
                uiZ.v64 &= ~lastBitMask;
            }
        } else if (softfloat_round_minMag != roundingMode) {
            if (is_sign(uiZ.v64) != (softfloat_round_max == roundingMode)) {
                uiZ.v64 = (uiZ.v64 | !!(uA.v0 != 0)) + roundBitsMask;
            }
        }

        uiZ.v64 &= ~roundBitsMask;
    }

    if (exact && (uiZ.v64 != uA.v64 || uiZ.v0 != uA.v0)) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return static_cast<float128_t>(uiZ);
}
