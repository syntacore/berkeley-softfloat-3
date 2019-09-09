
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

#include <utility>

namespace softfloat {
namespace internals {
namespace slow_int64 {

void
add_M_F128(uint32_t const* aWPtr,
           uint32_t const* bWPtr,
           uint32_t* const zWPtr,
           bool negateB)
{
    uint32_t uiA96 = aWPtr[index_word_hi(4)];
    int32_t expA = exp_F128_UI96(uiA96);
    uint32_t uiB96 = bWPtr[index_word_hi(4)];
    int32_t expB = exp_F128_UI96(uiB96);

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (try_propagate_NaN_M_F128(aWPtr, bWPtr, zWPtr)) {
            return;
        }

        uint32_t uiZ96 = uiA96;

        if (expB == 0x7FFF) {
            uiZ96 = uiB96 ^ pack_to_F128_UI96(negateB, 0, 0);

            if (0x7FFF == expA && uiZ96 != uiA96) {
                invalid_M_F128(zWPtr);
                return;
            }
        }

        zWPtr[index_word_hi(4)] = uiZ96;
        zWPtr[index_word(4, 2)] = 0;
        zWPtr[index_word(4, 1)] = 0;
        zWPtr[index_word(4, 0)] = 0;
        return;
    }

    bool signZ = is_sign(uiA96);
    bool const signB = is_sign(uiB96) != negateB;
    negateB = signZ != signB;

    if (static_cast<uint32_t>(uiA96 << 1) < static_cast<uint32_t>(uiB96 << 1)) {
        signZ = signB;
        expA = expB;
        expB = exp_F128_UI96(uiA96);
        std::swap(aWPtr, bWPtr);
        uiA96 = uiB96;
        uiB96 = bWPtr[index_word_hi(4)];
    }

    uint32_t sig96A = frac_F128_UI96(uiA96);
    uint32_t sig96B = frac_F128_UI96(uiB96);

    if (expA) {
        --expA;
        sig96A |= 0x00010000;

        if (expB) {
            --expB;
            sig96B |= 0x00010000;
        }
    }

    auto const addCarryMRoutinePtr = negateB ? add_compl_carry_M : add_carry_M;
    int32_t const expDiff = expA - expB;

    bool carry;
    uint32_t wordSigZ;
    uint32_t extSigZ[5];

    if (expDiff) {
        extSigZ[index_word_hi(5)] = sig96B;
        extSigZ[index_word(5, 3)] = bWPtr[index_word(4, 2)];
        extSigZ[index_word(5, 2)] = bWPtr[index_word(4, 1)];
        extSigZ[index_word(5, 1)] = bWPtr[index_word(4, 0)];
        extSigZ[index_word(5, 0)] = 0;
        shift_right_jam_M_160(extSigZ, static_cast<uint8_t>(expDiff), extSigZ);
        sig96B = extSigZ[index_word_hi(5)];
        carry = 0;

        if (negateB) {
            sig96B = ~sig96B;
            wordSigZ = extSigZ[index_word_lo(5)];
            extSigZ[index_word_lo(5)] = static_cast<uint32_t>(-static_cast<int32_t>(wordSigZ));
            carry = !wordSigZ;
        }

        carry = (*addCarryMRoutinePtr)(3,
                                       &aWPtr[index_multiword_lo(4, 3)],
                                       &extSigZ[index_multiword(5, 3, 1)],
                                       carry,
                                       &extSigZ[index_multiword(5, 3, 1)]);
        wordSigZ = sig96A + sig96B + !!(carry);
    } else {
        extSigZ[index_word_lo(5)] = 0;
        carry = (*addCarryMRoutinePtr)(3,
                                       &aWPtr[index_multiword_lo(4, 3)],
                                       &bWPtr[index_multiword_lo(4, 3)],
                                       negateB,
                                       &extSigZ[index_multiword(5, 3, 1)]);

        if (negateB) {
            wordSigZ = sig96A + ~sig96B + !!carry;

            if (wordSigZ & 0x80000000) {
                signZ = !signZ;
                carry = add_compl_carry_M_96(&bWPtr[index_multiword_lo(4, 3)],
                                             &aWPtr[index_multiword_lo(4, 3)],
                                             1,
                                             &extSigZ[index_multiword(5, 3, 1)]);
                wordSigZ = sig96B + ~sig96A + !!carry;
            } else {
                if (0 == wordSigZ && 0 == extSigZ[index_word(5, 3)] && 0 == (extSigZ[index_word(5, 2)] | extSigZ[index_word(5, 1)] | extSigZ[index_word(5, 0)])) {
                    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
                    signZ = softfloat_round_min == softfloat_roundingMode;
                    zWPtr[index_word_hi(4)] = pack_to_F128_UI96(signZ, 0, 0);
                    zWPtr[index_word(4, 2)] = 0;
                    zWPtr[index_word(4, 1)] = 0;
                    zWPtr[index_word(4, 0)] = 0;
                    return;
                }
            }
        } else {
            wordSigZ = sig96A + sig96B + !!carry;
        }
    }

    extSigZ[index_word_hi(5)] = wordSigZ;

    if (0x00010000 <= wordSigZ) {
        if (0x00020000 <= wordSigZ) {
            ++expA;
            short_shift_right_jam_M_160(extSigZ, 1, extSigZ);
        }

        round_pack_to_M_F128(signZ, expA, extSigZ, zWPtr);
    } else {
        norm_round_pack_to_M_F128(signZ, expA, extSigZ, zWPtr);
    }

}

}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
