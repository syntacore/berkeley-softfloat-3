
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015 The Regents of the University of
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

#include "internals.hpp"

#include "target.hpp"
#include "softfloat/functions.h"

namespace softfloat {
namespace internals {

extFloat80_t
softfloat_subMagsExtF80(uint16_t const uiA64,
                        uint64_t const uiA0,
                        uint16_t const uiB64,
                        uint64_t const uiB0,
                        bool const signZ)
{
    int32_t const expA = expExtF80UI64(uiA64);
    int32_t const expB = expExtF80UI64(uiB64);
    int32_t expDiff = expA - expB;

    if (0 < expDiff) {
        if (expA == 0x7FFF) {
            if (0 != (uiA0 & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
                uint128 const uiZ = softfloat_propagateNaNExtF80UI(uiA64, uiA0, uiB64, uiB0);
                extFloat80_t uZ;
                uZ.signExp = static_cast<uint16_t>(uiZ.v64);
                uZ.signif = uiZ.v0;
                return uZ;
            } else {
                extFloat80_t uZ;
                uZ.signExp = uiA64;
                uZ.signif = uiA0;
                return uZ;
            }
        } else {
            if (0 == expB) {
                --expDiff;

                if (0 == expDiff) {
                    uint128 const sig128_6 = softfloat_sub128(uiA0, 0, uiB0, 0);
                    return
                        softfloat_normRoundPackToExtF80(signZ, expA, sig128_6.v64, sig128_6.v0, extF80_roundingPrecision);
                }
            }

            uint128 const sig128_2 = softfloat_shiftRightJam128(uiB0, 0, static_cast<uint32_t>(expDiff));
            uint128 const sig128_1 = softfloat_sub128(uiA0, 0, sig128_2.v64, sig128_2.v0);
            return
                softfloat_normRoundPackToExtF80(signZ, expA, sig128_1.v64, sig128_1.v0, extF80_roundingPrecision);
        }
    } else if (expDiff < 0) {
        if (expB == 0x7FFF) {
            if (0 != (uiB0 & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
                uint128 const uiZ = softfloat_propagateNaNExtF80UI(uiA64, uiA0, uiB64, uiB0);
                extFloat80_t uZ;
                uZ.signExp = static_cast<uint16_t>(uiZ.v64);
                uZ.signif = uiZ.v0;
                return uZ;
            } else {
                extFloat80_t uZ;
                uZ.signExp = packToExtF80UI64(signZ ^ 1, 0x7FFF);
                uZ.signif = UINT64_C(0x8000000000000000);
                return uZ;
            }
        } else {
            if (0 == expA) {
                ++expDiff;

                if (0 == expDiff) {
                    uint128 const sig128_7 = softfloat_sub128(uiB0, 0, uiA0, 0);
                    return
                        softfloat_normRoundPackToExtF80(!signZ, expB, sig128_7.v64, sig128_7.v0, extF80_roundingPrecision);
                }
            }

            uint128 const sig128_4 = softfloat_shiftRightJam128(uiA0, 0, static_cast<uint32_t>(-expDiff));
            uint128 const sig128_3 = softfloat_sub128(uiB0, 0, sig128_4.v64, sig128_4.v0);
            return
                softfloat_normRoundPackToExtF80(!signZ, expB, sig128_3.v64, sig128_3.v0, extF80_roundingPrecision);
        }
    } else if (0x7FFF == expA) {
        if (0 != ((uiA0 | uiB0) & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
            uint128 const uiZ = softfloat_propagateNaNExtF80UI(uiA64, uiA0, uiB64, uiB0);
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(uiZ.v64);
            uZ.signif = uiZ.v0;
            return uZ;
        } else {
            softfloat_raiseFlags(softfloat_flag_invalid);
            extFloat80_t uZ;
            uZ.signExp = defaultNaNExtF80UI64;
            uZ.signif = defaultNaNExtF80UI0;
            return uZ;
        }
    } else {
        int32_t const expZ = 0 == expA ? 1 : expA;

        if (uiA0 > uiB0) {
            uint128 const sig128 = softfloat_sub128(uiA0, 0, uiB0, 0);
            return softfloat_normRoundPackToExtF80(signZ, expZ, sig128.v64, sig128.v0, extF80_roundingPrecision);
        } else if (uiA0 < uiB0) {
            uint128 const sig128 = softfloat_sub128(uiB0, 0, uiA0, 0);
            return softfloat_normRoundPackToExtF80(!signZ, expZ, sig128.v64, sig128.v0, extF80_roundingPrecision);
        } else {
            /* uiA0 == uiB0 */
            extFloat80_t uZ;
            uZ.signExp = packToExtF80UI64(softfloat_round_min == softfloat_roundingMode, 0);
            uZ.signif = 0;
            return uZ;
        }
    }
}

}  // namespace internals
}  // namespace softfloat
