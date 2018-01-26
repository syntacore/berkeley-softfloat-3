
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

#include "target.hpp"

namespace softfloat {
namespace internals {

float128_t
softfloat_subMagsF128(uint64_t uiA64,
                      uint64_t uiA0,
                      uint64_t uiB64,
                      uint64_t uiB0,
                      bool signZ)
{
    uint128 sigA;
    int32_t expB;
    uint128 sigB, sigZ;
    int32_t expDiff, expZ;
    uint128 uiZ;

    int32_t const expA = expF128UI64(uiA64);
    sigA.v64 = fracF128UI64(uiA64);
    sigA.v0 = uiA0;
    expB = expF128UI64(uiB64);
    sigB.v64 = fracF128UI64(uiB64);
    sigB.v0 = uiB0;
    sigA = softfloat_shortShiftLeft128(sigA.v64, sigA.v0, 4);
    sigB = softfloat_shortShiftLeft128(sigB.v64, sigB.v0, 4);
    expDiff = expA - expB;

    if (0 < expDiff) {
        goto expABigger;
    }

    if (expDiff < 0) {
        goto expBBigger;
    }

    if (expA == 0x7FFF) {
        if (sigA.v64 | sigA.v0 | sigB.v64 | sigB.v0) {
            goto propagateNaN;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        uiZ.v64 = defaultNaNF128UI64;
        uiZ.v0 = defaultNaNF128UI0;
        goto uiZ;
    }

    expZ = expA;

    if (!expZ) {
        expZ = 1;
    }

    if (sigB.v64 < sigA.v64) {
        goto aBigger;
    } else if (sigA.v64 < sigB.v64) {
        goto bBigger;
    } else if (sigB.v0 < sigA.v0) {
        goto aBigger;
    } else if (sigA.v0 < sigB.v0) {
        goto bBigger;
    } else {
        uiZ.v64 = packToF128UI64(softfloat_roundingMode == softfloat_round_min, 0, 0);
        uiZ.v0 = 0;
        goto uiZ;
    }

expBBigger:

    if (expB == 0x7FFF) {
        if (sigB.v64 | sigB.v0) {
            goto propagateNaN;
        }

        uiZ.v64 = packToF128UI64(signZ ^ 1, 0x7FFF, 0);
        uiZ.v0 = 0;
        goto uiZ;
    }

    if (expA) {
        sigA.v64 |= UINT64_C(0x0010000000000000);
    } else {
        ++expDiff;

        if (!expDiff) {
            goto newlyAlignedBBigger;
        }
    }

    sigA = softfloat_shiftRightJam128(sigA.v64, sigA.v0, static_cast<uint32_t>(-expDiff));
newlyAlignedBBigger:
    expZ = expB;
    sigB.v64 |= UINT64_C(0x0010000000000000);
bBigger:
    signZ = !signZ;
    sigZ = softfloat_sub128(sigB.v64, sigB.v0, sigA.v64, sigA.v0);
    goto normRoundPack;
expABigger:

    if (expA == 0x7FFF) {
        if (sigA.v64 | sigA.v0) {
            goto propagateNaN;
        }

        uiZ.v64 = uiA64;
        uiZ.v0 = uiA0;
        goto uiZ;
    }

    if (expB) {
        sigB.v64 |= UINT64_C(0x0010000000000000);
    } else {
        --expDiff;

        if (!expDiff) {
            goto newlyAlignedABigger;
        }
    }

    sigB = softfloat_shiftRightJam128(sigB.v64, sigB.v0, static_cast<uint32_t>(expDiff));
newlyAlignedABigger:
    expZ = expA;
    sigA.v64 |= UINT64_C(0x0010000000000000);
aBigger:
    sigZ = softfloat_sub128(sigA.v64, sigA.v0, sigB.v64, sigB.v0);
normRoundPack:
    return softfloat_normRoundPackToF128(signZ, expZ - 5, sigZ.v64, sigZ.v0);
propagateNaN:
    uiZ = softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0);
uiZ:
    return u_as_f_128(uiZ);
}

}  // namespace internals
}  // namespace softfloat
