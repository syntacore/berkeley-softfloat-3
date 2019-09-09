
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

void
norm_round_pack_to_M_F128(bool const sign,
                          int32_t exp,
                          uint32_t* const extSigPtr,
                          uint32_t* const zWPtr)
{
    auto ptr = extSigPtr + index_word_hi(5);
    int16_t shiftDist = 0;
    uint32_t wordSig;

    for (;;) {
        wordSig = *ptr;

        if (wordSig) {
            break;
        }

        shiftDist += 32;

        if (160 <= shiftDist) {
            zWPtr[index_word_hi(4)] = pack_to_F128_UI96(sign, 0, 0);
            zWPtr[index_word(4, 2)] = 0;
            zWPtr[index_word(4, 1)] = 0;
            zWPtr[index_word(4, 0)] = 0;
            return;
        }

        ptr -= wordIncr;
    }

    int16_t const shiftDist_1 = shiftDist + count_leading_zeros(wordSig) - 15;

    if (0 != shiftDist_1) {
        exp -= shiftDist_1;
        /**
        @bug modify input
        */
        shift_left_M_160(extSigPtr, static_cast<uint32_t>(shiftDist_1), extSigPtr);
    }

    round_pack_to_M_F128(sign, exp, extSigPtr, zWPtr);
}

}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
