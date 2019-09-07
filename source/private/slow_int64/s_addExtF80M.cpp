
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
#include <algorithm>

namespace softfloat {
namespace internals {
namespace slow_int64 {

/**
@bug use extFloat80_t
*/
void
softfloat_addExtF80M(extFloat80M const* aSPtr,
                     extFloat80M const* bSPtr,
                     extFloat80M* const zSPtr,
                     bool negateB)
{
    typedef void(*round_packer_ptr)(bool, int32_t, uint32_t*, uint8_t const&, extFloat80M*);

    uint16_t const uiA64 = aSPtr->signExp;
    int32_t expA = expExtF80UI64(uiA64);
    uint16_t const uiB64 = bSPtr->signExp;
    int32_t expB = expExtF80UI64(uiB64);

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (softfloat_tryPropagateNaNExtF80M(aSPtr, bSPtr, zSPtr)) {
            return;
        }

        uint16_t uiZ64 = uiA64;

        if (0x7FFF == expB) {
            uiZ64 = static_cast<uint16_t>(uiB64 ^ packToExtF80UI64(negateB, 0));

            if (0x7FFF == expA && uiZ64 != uiA64) {
                softfloat_invalidExtF80M(zSPtr);
                return;
            }
        }

        zSPtr->signExp = uiZ64;
        zSPtr->signif = static_cast<uint64_t>(INT64_MIN);
        return;
    }

    bool signZ = is_sign(uiA64);
    bool const signB = is_sign(uiB64) != negateB;
    negateB = signZ != signB;

    if (expA < expB) {
        signZ = signB;
        std::swap(expA, expB);
        std::swap(aSPtr, bSPtr);
    }

    if (0 == expB) {
        expB = 1;

        if (0 == expA) {
            expA = 1;
        }
    }

    uint64_t sigZ = aSPtr->signif;
    uint64_t sigB = bSPtr->signif;
    round_packer_ptr roundPackRoutinePtr = softfloat_roundPackMToExtF80M;
    int32_t const expDiff = expA - expB;

    if (0 != expDiff) {
        uint32_t extSigX[3];
        extSigX[indexWord(3, 2)] = sigB >> 32;
        extSigX[indexWord(3, 1)] = static_cast<uint32_t>(sigB);
        extSigX[indexWord(3, 0)] = 0;
        softfloat_shiftRightJam96M(extSigX, static_cast<uint8_t>(expDiff), extSigX);
        sigB =
            static_cast<uint64_t>(extSigX[indexWord(3, 2)]) << 32 |
            extSigX[indexWord(3, 1)];

        if (!negateB) {
            sigZ += sigB;

            if (0 == (sigZ & UINT64_C(0x8000000000000000))) {
                auto const sigZExtra_1 = (static_cast<uint32_t>(sigZ) << 31) | !!(0 != extSigX[indexWordLo(3)]);
                ++expA;
                sigZ = UINT64_C(0x8000000000000000) | sigZ >> 1;
                extSigX[indexWord(3, 0)] = sigZExtra_1;
            }

            extSigX[indexWord(3, 2)] = sigZ >> 32;
            extSigX[indexWord(3, 1)] = static_cast<uint32_t>(sigZ);
            (*roundPackRoutinePtr)(signZ, expA, extSigX, extF80_roundingPrecision, zSPtr);
            return;
        }

        sigZ -= sigB;
        uint32_t sigZExtra = extSigX[indexWordLo(3)];

        if (0 != sigZExtra) {
            --sigZ;
            sigZExtra = static_cast<uint16_t>(-static_cast<int32_t>(sigZExtra));
        }

        if (0 == (sigZ & UINT64_C(0x8000000000000000))) {
            if (0 != (sigZ & UINT64_C(0x4000000000000000))) {
                --expA;
                sigZ = sigZ << 1 | sigZExtra >> 31;
                sigZExtra <<= 1;
            } else {
                roundPackRoutinePtr = softfloat_normRoundPackMToExtF80M;
            }
        }

        {
            uint32_t extSigX[3];
            extSigX[indexWord(3, 0)] = sigZExtra;
            extSigX[indexWord(3, 2)] = sigZ >> 32;
            extSigX[indexWord(3, 1)] = static_cast<uint32_t>(sigZ);
            (*roundPackRoutinePtr)(signZ, expA, extSigX, extF80_roundingPrecision, zSPtr);
        }

        return;
    }

    assert(0 == expDiff);
    uint32_t sigZExtra = 0;

    if (negateB) {
        if (sigZ < sigB) {
            signZ = !signZ;
            sigZ = sigB - sigZ;
        } else {
            sigZ -= sigB;

            if (0 == sigZ) {
                softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
                zSPtr->signExp = packToExtF80UI64(softfloat_round_min == softfloat_roundingMode, 0);
                zSPtr->signif = 0;
                return;
            }
        }

        roundPackRoutinePtr = softfloat_normRoundPackMToExtF80M;
    } else {
        sigZ += sigB;

        if (sigZ < sigB) {
            sigZExtra = static_cast<uint32_t>(sigZ) << 31;
            ++expA;
            sigZ = UINT64_C(0x8000000000000000) | sigZ >> 1;
        } else if (0 == (sigZ & UINT64_C(0x8000000000000000))) {
            roundPackRoutinePtr = softfloat_normRoundPackMToExtF80M;
        }
    }

    {
        uint32_t extSigX[3];
        extSigX[indexWord(3, 0)] = sigZExtra;
        extSigX[indexWord(3, 2)] = sigZ >> 32;
        extSigX[indexWord(3, 1)] = static_cast<uint32_t>(sigZ);
        (*roundPackRoutinePtr)(signZ, expA, extSigX, extF80_roundingPrecision, zSPtr);
    }
}

}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
