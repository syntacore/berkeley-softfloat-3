
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

#include "model.hpp"

#ifndef SOFTFLOAT_FAST_INT64
#error Fast int64_t operations only
#endif

namespace softfloat {
namespace internals {

float128_t
softfloat_subMagsF128(uint64_t const uiA64,
                      uint64_t const uiA0,
                      uint64_t const uiB64,
                      uint64_t const uiB0,
                      bool const signZ)
{
    int32_t const expA = expF128UI64(uiA64);
    int32_t const expB = expF128UI64(uiB64);
    uint128 const sigA = softfloat_shortShiftLeft128(fracF128UI64(uiA64), uiA0, 4);
    uint128 const sigB = softfloat_shortShiftLeft128(fracF128UI64(uiB64), uiB0, 4);
    int32_t const expDiff = expA - expB;

    if (!(0 < expDiff)) {
        if (!(expDiff < 0)) {

            if (0x7FFF == expA) {
                if (sigA.v64 | sigA.v0 | sigB.v64 | sigB.v0) {
                    return float128_t(softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0));
                }

                softfloat_raiseFlags(softfloat_flag_invalid);
                return float128_t(uint128{defaultNaNF128UI64, defaultNaNF128UI0});
            }

            int32_t const expZ = 0 == expA ? 1 : expA;

            if (sigB.v64 < sigA.v64) {
                uint128 const sigZ = softfloat_sub128(sigA.v64, sigA.v0, sigB.v64, sigB.v0);
                return softfloat_normRoundPackToF128(signZ, expZ - 5, sigZ.v64, sigZ.v0);
            }

            if (sigA.v64 < sigB.v64) {
                uint128 const sigZ = softfloat_sub128(sigB.v64, sigB.v0, sigA.v64, sigA.v0);
                return softfloat_normRoundPackToF128(!signZ, expZ - 5, sigZ.v64, sigZ.v0);
            }

            if (sigB.v0 < sigA.v0) {
                uint128 const sigZ = softfloat_sub128(sigA.v64, sigA.v0, sigB.v64, sigB.v0);
                return softfloat_normRoundPackToF128(signZ, expZ - 5, sigZ.v64, sigZ.v0);
            }

            if (sigA.v0 < sigB.v0) {
                uint128 const sigZ = softfloat_sub128(sigB.v64, sigB.v0, sigA.v64, sigA.v0);
                return softfloat_normRoundPackToF128(!signZ, expZ - 5, sigZ.v64, sigZ.v0);
            }

            softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
            return float128_t(uint128{packToF128UI64(softfloat_roundingMode == softfloat_round_min, 0, 0), 0});
        }

        if (0x7FFF == expB) {
            if (0 != (sigB.v64 | sigB.v0)) {
                return float128_t(softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0));
            }

            return float128_t(uint128{packToF128UI64(!signZ, 0x7FFF, 0), 0});
        }

        if (0 != expA) {
            auto const sigA_1 = softfloat_shiftRightJam128(sigA.v64 | UINT64_C(0x0010000000000000), sigA.v0, static_cast<uint32_t>(-expDiff));
            uint128 const sigZ = softfloat_sub128(sigB.v64 | UINT64_C(0x0010000000000000), sigB.v0, sigA_1.v64, sigA_1.v0);
            return softfloat_normRoundPackToF128(!signZ, expB - 5, sigZ.v64, sigZ.v0);
        }

        auto const expDiff_1 = expDiff + 1;

        if (0 == expDiff_1) {
            uint128 const sigZ = softfloat_sub128(sigB.v64 | UINT64_C(0x0010000000000000), sigB.v0, sigA.v64, sigA.v0);
            return softfloat_normRoundPackToF128(!signZ, expB - 5, sigZ.v64, sigZ.v0);
        }

        auto const sigA_1 = softfloat_shiftRightJam128(sigA.v64, sigA.v0, static_cast<uint32_t>(-expDiff_1));
        uint128 const sigZ = softfloat_sub128(sigB.v64 | UINT64_C(0x0010000000000000), sigB.v0, sigA_1.v64, sigA_1.v0);
        return softfloat_normRoundPackToF128(!signZ, expB - 5, sigZ.v64, sigZ.v0);
    }

    if (0x7FFF == expA) {
        if (0 != (sigA.v64 | sigA.v0)) {
            return float128_t(softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0));
        }

        return float128_t(uint128{uiA64, uiA0});
    }

    if (0 != expB) {
        auto const sigB_1 = softfloat_shiftRightJam128(sigB.v64 | UINT64_C(0x0010000000000000), sigB.v0, static_cast<uint32_t>(expDiff));
        uint128 const sigZ = softfloat_sub128(sigA.v64 | UINT64_C(0x0010000000000000), sigA.v0, sigB_1.v64, sigB_1.v0);
        return softfloat_normRoundPackToF128(signZ, expA - 5, sigZ.v64, sigZ.v0);
    }

    auto const expDiff_1 = expDiff - 1;

    if (0 == expDiff_1) {
        uint128 const sigZ = softfloat_sub128(sigA.v64 | UINT64_C(0x0010000000000000), sigA.v0, sigB.v64, sigB.v0);
        return softfloat_normRoundPackToF128(signZ, expA - 5, sigZ.v64, sigZ.v0);
    }

    auto const sigB_1 = softfloat_shiftRightJam128(sigB.v64, sigB.v0, static_cast<uint32_t>(expDiff_1));
    uint128 const sigZ = softfloat_sub128(sigA.v64 | UINT64_C(0x0010000000000000), sigA.v0, sigB_1.v64, sigB_1.v0);
    return softfloat_normRoundPackToF128(signZ, expA - 5, sigZ.v64, sigZ.v0);
}

}  // namespace internals
}  // namespace softfloat
