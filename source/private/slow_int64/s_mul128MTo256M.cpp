
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
mul_M_128_to_256(uint32_t const* aPtr,
                 uint32_t const* bPtr,
                 uint32_t* zPtr)
{
    bPtr += index_word_lo(4);
    uint32_t* const lastZPtr = zPtr + index_multiword_hi(8, 5);
    zPtr += index_multiword_lo(8, 5);
    uint32_t wordB = *bPtr;
    uint64_t dwordProd = static_cast<uint64_t>(aPtr[index_word(4, 0)]) * wordB;
    zPtr[index_word(5, 0)] = static_cast<uint32_t>(dwordProd);
    dwordProd = static_cast<uint64_t>(aPtr[index_word(4, 1)]) * wordB + (dwordProd >> 32);
    zPtr[index_word(5, 1)] = static_cast<uint32_t>(dwordProd);
    dwordProd = static_cast<uint64_t>(aPtr[index_word(4, 2)]) * wordB + (dwordProd >> 32);
    zPtr[index_word(5, 2)] = static_cast<uint32_t>(dwordProd);
    dwordProd = static_cast<uint64_t>(aPtr[index_word(4, 3)]) * wordB + (dwordProd >> 32);
    zPtr[index_word(5, 3)] = static_cast<uint32_t>(dwordProd);
    zPtr[index_word(5, 4)] = dwordProd >> 32;

    do {
        bPtr += wordIncr;
        zPtr += wordIncr;
        wordB = *bPtr;
        dwordProd = static_cast<uint64_t>(aPtr[index_word(4, 0)]) * wordB;
        uint32_t wordZ = zPtr[index_word(5, 0)] + static_cast<uint32_t>(dwordProd);
        zPtr[index_word(5, 0)] = wordZ;
        bool carry = wordZ < static_cast<uint32_t>(dwordProd);
        dwordProd =
            static_cast<uint64_t>(aPtr[index_word(4, 1)]) * wordB + (dwordProd >> 32);
        wordZ = zPtr[index_word(5, 1)] + static_cast<uint32_t>(dwordProd) + !!carry;
        zPtr[index_word(5, 1)] = wordZ;

        if (wordZ != static_cast<uint32_t>(dwordProd)) {
            carry = wordZ < static_cast<uint32_t>(dwordProd);
        }

        dwordProd =
            static_cast<uint64_t>(aPtr[index_word(4, 2)]) * wordB + (dwordProd >> 32);
        wordZ = zPtr[index_word(5, 2)] + static_cast<uint32_t>(dwordProd) + !!carry;
        zPtr[index_word(5, 2)] = wordZ;

        if (wordZ != static_cast<uint32_t>(dwordProd)) {
            carry = wordZ < static_cast<uint32_t>(dwordProd);
        }

        dwordProd =
            static_cast<uint64_t>(aPtr[index_word(4, 3)]) * wordB + (dwordProd >> 32);
        wordZ = zPtr[index_word(5, 3)] + static_cast<uint32_t>(dwordProd) + !!carry;
        zPtr[index_word(5, 3)] = wordZ;

        if (wordZ != static_cast<uint32_t>(dwordProd)) {
            carry = wordZ < static_cast<uint32_t>(dwordProd);
        }

        zPtr[index_word(5, 4)] = (dwordProd >> 32) + !!carry;
    } while (zPtr != lastZPtr);
}

}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
