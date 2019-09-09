
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
shift_left_M(size_t size_words,
                     uint32_t const* aPtr,
                     uint32_t dist,
                     uint32_t* zPtr)
{
    uint32_t* destPtr;
    uint32_t wordDist = dist >> 5;

    if (wordDist < size_words) {
        aPtr += index_multiword_lo_but(size_words, wordDist);
        uint8_t const innerDist = dist & 31;

        if (innerDist) {
            short_shift_left_M(static_cast<uint8_t>(size_words - wordDist),
                                      aPtr,
                                      innerDist,
                                      zPtr + index_multiword_hi_but(size_words, wordDist));

            if (0 == wordDist) {
                return;
            }
        } else {
            aPtr += index_word_hi(size_words - wordDist);
            destPtr = zPtr + index_word_hi(size_words);

            for (size_t i = size_words - wordDist; i; --i) {
                *destPtr = *aPtr;
                aPtr -= wordIncr;
                destPtr -= wordIncr;
            }
        }

        zPtr += index_multiword_lo(size_words, wordDist);
    } else {
        wordDist = size_words;
    }

    do {
        *zPtr++ = 0;
        --wordDist;
    } while (wordDist);

}

}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
