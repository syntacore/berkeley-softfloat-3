
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

float128_t
softfloat_addMagsF128(uint64_t const uiA64,
                      uint64_t const uiA0,
                      uint64_t const uiB64,
                      uint64_t const uiB0,
                      bool const signZ)
{
    int32_t const expA = expF128UI64(uiA64);
    uint128 sigA{fracF128UI64(uiA64), uiA0};
    int32_t const expB = expF128UI64(uiB64);
    uint128 sigB{fracF128UI64(uiB64), uiB0};
    int32_t expDiff = expA - expB;

    if (0 == expDiff) {
        if (0x7FFF == expA) {
            if (0 != (sigA.v64 | sigA.v0 | sigB.v64 | sigB.v0)) {
                return static_cast<float128_t>(softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0));
            }

            return static_cast<float128_t>(uint128{uiA64, uiA0});
        }

        uint128 const sigZ = softfloat_add128(sigA.v64, sigA.v0, sigB.v64, sigB.v0);

        if (0 == expA) {
            return static_cast<float128_t>(uint128{packToF128UI64(signZ, 0, sigZ.v64), sigZ.v0});
        }

        uint128_extra const sig128Extra = softfloat_shortShiftRightJam128Extra(sigZ.v64 | UINT64_C(0x0002000000000000), sigZ.v0, 0, 1);
        return softfloat_roundPackToF128(signZ, expA, sig128Extra.v.v64, sig128Extra.v.v0, sig128Extra.extra);
    }

    if (expDiff < 0) {
        if (0x7FFF == expB) {
            if (sigB.v64 | sigB.v0) {
                return static_cast<float128_t>(softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0));
            }

            return static_cast<float128_t>(uint128{packToF128UI64(signZ, 0x7FFF, 0), 0});
        }

        if (0 != expA) {
            sigA.v64 |= UINT64_C(0x0001000000000000);
        } else if (0 == expDiff + 1) {
            uint128 const sigZ = softfloat_add128(sigA.v64 | UINT64_C(0x0001000000000000), sigA.v0, sigB.v64, sigB.v0);

            if (sigZ.v64 < UINT64_C(0x0002000000000000)) {
                return softfloat_roundPackToF128(signZ, expB - 1, sigZ.v64, sigZ.v0, 0);
            }

            uint128_extra const sig128Extra = softfloat_shortShiftRightJam128Extra(sigZ.v64, sigZ.v0, 0, 1);
            return softfloat_roundPackToF128(signZ, expB, sig128Extra.v.v64, sig128Extra.v.v0, sig128Extra.extra);
        }

        uint128_extra const sig128Extra = softfloat_shiftRightJam128Extra(sigA.v64, sigA.v0, 0, static_cast<uint32_t>(-(expDiff + 1)));
        uint128 const sigZ = softfloat_add128(sig128Extra.v.v64 | UINT64_C(0x0001000000000000), sig128Extra.v.v0, sigB.v64, sigB.v0);

        if (sigZ.v64 < UINT64_C(0x0002000000000000)) {
            return softfloat_roundPackToF128(signZ, expB - 1, sigZ.v64, sigZ.v0, sig128Extra.extra);
        }

        uint128_extra const sig128Extra_1 = softfloat_shortShiftRightJam128Extra(sigZ.v64, sigZ.v0, sig128Extra.extra, 1);
        return softfloat_roundPackToF128(signZ, expB, sig128Extra_1.v.v64, sig128Extra_1.v.v0, sig128Extra_1.extra);
    }

    if (0x7FFF == expA) {
        if (0 != (sigA.v64 | sigA.v0)) {
            return static_cast<float128_t>(softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0));
        }

        return static_cast<float128_t>(uint128{uiA64, uiA0});
    }

    if (0 != expB) {
        sigB.v64 |= UINT64_C(0x0001000000000000);
    } else if (1 == expDiff) {
        uint128 const sigZ = softfloat_add128(sigA.v64 | UINT64_C(0x0001000000000000), sigA.v0, sigB.v64, sigB.v0);

        if (sigZ.v64 < UINT64_C(0x0002000000000000)) {
            return softfloat_roundPackToF128(signZ, expA - 1, sigZ.v64, sigZ.v0, 0);
        }

        uint128_extra const sig128Extra = softfloat_shortShiftRightJam128Extra(sigZ.v64, sigZ.v0, 0, 1);
        return softfloat_roundPackToF128(signZ, expA, sig128Extra.v.v64, sig128Extra.v.v0, sig128Extra.extra);
    }

    uint128_extra const sig128Extra = softfloat_shiftRightJam128Extra(sigB.v64, sigB.v0, 0, static_cast<uint32_t>(expDiff - 1));
    uint128 const sigZ = softfloat_add128(sigA.v64 | UINT64_C(0x0001000000000000), sigA.v0, sig128Extra.v.v64, sig128Extra.v.v0);

    if (sigZ.v64 < UINT64_C(0x0002000000000000)) {
        return softfloat_roundPackToF128(signZ, expA - 1, sigZ.v64, sigZ.v0, sig128Extra.extra);
    }

    uint128_extra const sig128Extra_1 = softfloat_shortShiftRightJam128Extra(sigZ.v64, sigZ.v0, sig128Extra.extra, 1);
    return softfloat_roundPackToF128(signZ, expA, sig128Extra_1.v.v64, sig128Extra_1.v.v0, sig128Extra_1.extra);
}

}  // namespace internals
}  // namespace softfloat
