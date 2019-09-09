
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

namespace softfloat {
namespace internals {
namespace fast_int64 {

extFloat80_t
add_magnitudes(uint16_t const uiA64,
               uint64_t const uiA0,
               uint16_t const uiB64,
               uint64_t const uiB0,
               bool signZ)
{
    int32_t const expA = exp_extF80_UI64(uiA64);
    int32_t const expB = exp_extF80_UI64(uiB64);
    int32_t expDiff = expA - expB;

    if (0 == expDiff) {
        if (0x7FFF == expA) {
            extFloat80_t uZ;

            if (0 != (UINT64_C(0x7FFFFFFFFFFFFFFF) & (uiA0 | uiB0))) {
                uint128 const uiZ = propagate_NaN(uiA64, uiA0, uiB64, uiB0);
                uZ.signExp = static_cast<uint16_t>(uiZ.v64);
                uZ.signif = uiZ.v0;
            } else {
                uZ.signExp = uiA64;
                uZ.signif = uiA0;
            }

            return uZ;
        }

        auto const sigZ_1 = uiA0 + uiB0;

        if (0 == expA) {
            exp32_sig64 const normExpSig = norm_subnormal_extF80Sig(sigZ_1);
            return round_pack_to_extF80(signZ, normExpSig.exp + 1, normExpSig.sig, 0, extF80_rounding_precision);
        }

        uint64_extra const sig64Extra = shortShiftRightJam64Extra(sigZ_1, 0, 1);
        return round_pack_to_extF80(signZ, expA + 1, UINT64_C(0x8000000000000000) | sig64Extra.v, sig64Extra.extra, extF80_rounding_precision);
    }

    if (expDiff < 0) {
        if (0x7FFF == expB) {
            extFloat80_t uZ;

            if (uiB0 & UINT64_C(0x7FFFFFFFFFFFFFFF)) {
                uint128 const uiZ = propagate_NaN(uiA64, uiA0, uiB64, uiB0);
                uZ.signExp = static_cast<uint16_t>(uiZ.v64);
                uZ.signif = uiZ.v0;
            } else {
                uZ.signExp = pack_to_extF80_UI64(signZ, 0x7FFF);
                uZ.signif = uiB0;
            }

            return uZ;
        }

        auto const expZ = expB;

        if (0 == expA) {
            ++expDiff;

            if (0 == expDiff) {
                auto const sigZ_1 = uiA0 + uiB0;

                if (0 != (UINT64_C(0x8000000000000000) & sigZ_1)) {
                    return round_pack_to_extF80(signZ, expZ, sigZ_1, 0, extF80_rounding_precision);
                }

                uint64_extra const sig64Extra = shortShiftRightJam64Extra(sigZ_1, 0, 1);
                return round_pack_to_extF80(signZ, expZ + 1, UINT64_C(0x8000000000000000) | sig64Extra.v, sig64Extra.extra, extF80_rounding_precision);
            }
        }

        uint64_extra const sig64Extra = shift_right_jam_64Extra(uiA0, 0, static_cast<uint32_t>(-expDiff));
        auto const sigZ_1 = sig64Extra.v + uiB0;

        if (0 != (sigZ_1 & UINT64_C(0x8000000000000000))) {
            return round_pack_to_extF80(signZ, expZ, sigZ_1, sig64Extra.extra, extF80_rounding_precision);
        }

        auto const sig64Extra_1 = shortShiftRightJam64Extra(sigZ_1, sig64Extra.extra, 1);
        return round_pack_to_extF80(signZ, expZ + 1, sig64Extra_1.v | UINT64_C(0x8000000000000000), sig64Extra_1.extra, extF80_rounding_precision);
    }

    if (0x7FFF == expA) {
        if (0 != (UINT64_C(0x7FFFFFFFFFFFFFFF) & uiA0)) {
            uint128 const uiZ = propagate_NaN(uiA64, uiA0, uiB64, uiB0);
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
    }

    auto const expZ = expA;

    if (0 == expB) {
        --expDiff;

        if (0 == expDiff) {
            auto const sigZ_1 = uiA0 + uiB0;

            if (0 != (UINT64_C(0x8000000000000000) & sigZ_1)) {
                return round_pack_to_extF80(signZ, expZ, sigZ_1, 0, extF80_rounding_precision);
            }

            uint64_extra const sig64Extra = shortShiftRightJam64Extra(sigZ_1, 0, 1);
            return round_pack_to_extF80(signZ, expZ + 1, UINT64_C(0x8000000000000000) | sig64Extra.v, sig64Extra.extra, extF80_rounding_precision);
        }
    }

    uint64_extra const sig64Extra = shift_right_jam_64Extra(uiB0, 0u, static_cast<uint32_t>(expDiff));
    auto const sigZ_1 = uiA0 + sig64Extra.v;

    if (0 != (sigZ_1 & UINT64_C(0x8000000000000000))) {
        return round_pack_to_extF80(signZ, expZ, sigZ_1, sig64Extra.extra, extF80_rounding_precision);
    }

    auto const sig64Extra_1 = shortShiftRightJam64Extra(sigZ_1, sig64Extra.extra, 1);
    return round_pack_to_extF80(signZ, expZ + 1, sig64Extra_1.v | UINT64_C(0x8000000000000000), sig64Extra_1.extra, extF80_rounding_precision);
}

}  // namespace fast_int64
}  // namespace internals
}  // namespace softfloat
