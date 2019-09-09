
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

#ifdef BIG_ENDIAN
#define INIT_UINTM4( v3, v2, v1, v0 ) { v3, v2, v1, v0 }
#else
#define INIT_UINTM4( v3, v2, v1, v0 ) { v0, v1, v2, v3 }
#endif  // BIG_ENDIAN

namespace softfloat {
namespace internals {
namespace slow_int64 {

void
round_pack_to_M_F128(bool const sign,
                            int32_t exp,
                            uint32_t* const extSigPtr,
                            uint32_t* const zWPtr)
{
    static const uint32_t maxSig[4] = INIT_UINTM4(0x0001FFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);

    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    bool const roundNearEven = softfloat_round_near_even == softfloat_roundingMode;
    uint32_t sigExtra = extSigPtr[index_word_lo(5)];
    bool doIncrement = (0x80000000 <= sigExtra);

    if (!roundNearEven && softfloat_round_near_maxMag != softfloat_roundingMode) {
        doIncrement =
            (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode &&
            sigExtra;
    }

    if (0x7FFD <= static_cast<uint32_t>(exp)) {
        if (exp < 0) {
            bool const isTiny =
                softfloat_tininess_beforeRounding == detect_tininess ||
                exp < -1 ||
                !doIncrement ||
                compare_M_128(extSigPtr + index_multiword_hi(5, 4), maxSig) < 0;

            /**
            @bug modify input
            */
            shift_right_jam_M_160(extSigPtr, static_cast<uint8_t>(-exp), extSigPtr);
            exp = 0;
            sigExtra = extSigPtr[index_word_lo(5)];

            if (isTiny && sigExtra) {
                softfloat_raiseFlags(softfloat_flag_underflow);
            }

            doIncrement = 0x80000000 <= sigExtra;

            if (!roundNearEven && softfloat_round_near_maxMag != softfloat_roundingMode) {
                doIncrement =
                    ((sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode) &&
                    sigExtra;
            }
        } else if (0x7FFD < exp || (0x7FFD == exp && doIncrement && 0 == compare_M_128(extSigPtr + index_multiword_hi(5, 4), maxSig))) {
            softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);

            if (roundNearEven || softfloat_round_near_maxMag == softfloat_roundingMode || (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode) {
                zWPtr[index_word_hi(4)] = pack_to_F128_UI96(sign, 0x7FFF, 0);
                zWPtr[index_word(4, 2)] = 0;
                zWPtr[index_word(4, 1)] = 0;
                zWPtr[index_word(4, 0)] = 0;
                return;
            }

            zWPtr[index_word_hi(4)] = pack_to_F128_UI96(sign, 0x7FFE, 0x0000FFFF);
            zWPtr[index_word(4, 2)] = 0xFFFFFFFF;
            zWPtr[index_word(4, 1)] = 0xFFFFFFFF;
            zWPtr[index_word(4, 0)] = 0xFFFFFFFF;
            return;
        }
    }

    if (sigExtra) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    uint32_t uj = extSigPtr[index_word(5, 1)];

    if (doIncrement) {
        uint32_t ui;
        ++uj;

        if (uj) {
            if (!(sigExtra & 0x7FFFFFFF) && roundNearEven) {
                uj &= ~1;
            }

            zWPtr[index_word(4, 2)] = extSigPtr[index_word(5, 3)];
            zWPtr[index_word(4, 1)] = extSigPtr[index_word(5, 2)];
            zWPtr[index_word(4, 0)] = uj;
            ui = extSigPtr[index_word_hi(5)];
        } else {
            zWPtr[index_word(4, 0)] = uj;
            ui = extSigPtr[index_word(5, 2)] + 1;
            zWPtr[index_word(4, 1)] = ui;
            uj = extSigPtr[index_word(5, 3)];

            if (ui) {
                zWPtr[index_word(4, 2)] = uj;
                ui = extSigPtr[index_word_hi(5)];
            } else {
                ++uj;
                zWPtr[index_word(4, 2)] = uj;
                ui = extSigPtr[index_word_hi(5)];

                if (!uj) {
                    ++ui;
                }
            }
        }

        zWPtr[index_word_hi(4)] = pack_to_F128_UI96(sign, static_cast<unsigned>(exp), ui);
    } else {
        zWPtr[index_word(4, 0)] = uj;
        uint32_t ui = extSigPtr[index_word(5, 2)];
        zWPtr[index_word(4, 1)] = ui;
        uj |= ui;
        ui = extSigPtr[index_word(5, 3)];
        zWPtr[index_word(4, 2)] = ui;
        uj |= ui;
        ui = extSigPtr[index_word_hi(5)];
        uj |= ui;

        if (0 == uj) {
            exp = 0;
        }

        zWPtr[index_word_hi(4)] = pack_to_F128_UI96(sign, static_cast<unsigned>(exp), ui);
    }
}

}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
