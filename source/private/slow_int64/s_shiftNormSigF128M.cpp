
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
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

namespace softfloat {
namespace internals {
namespace slow_int64 {

int
shift_norm_sig_M_F128(const uint32_t* wPtr,
                     uint8_t shiftDist,
                     uint32_t* sigPtr)
{
    uint32_t wordSig = wPtr[index_word_hi(4)];
    int32_t exp = exp_F128_UI96(wordSig);

    if (exp) {
        short_shift_left_M_128(wPtr, shiftDist, sigPtr);
        uint32_t leadingBit = UINT32_C(0x00010000) << shiftDist;
        sigPtr[index_word_hi(4)] =
            (sigPtr[index_word_hi(4)] & (leadingBit - 1)) | leadingBit;
    } else {
        exp = 16;
        wordSig &= 0x7FFFFFFF;

        if (!wordSig) {
            exp = -16;
            wordSig = wPtr[index_word(4, 2)];

            if (!wordSig) {
                exp = -48;
                wordSig = wPtr[index_word(4, 1)];

                if (!wordSig) {
                    wordSig = wPtr[index_word(4, 0)];

                    if (!wordSig) {
                        return -128;
                    }

                    exp = -80;
                }
            }
        }

        exp -= count_leading_zeros(wordSig);
        shift_left_M_128(wPtr, static_cast<uint32_t>(1 - exp + shiftDist), sigPtr);
    }

    return exp;

}

}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
