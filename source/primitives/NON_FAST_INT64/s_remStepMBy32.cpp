
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

void
softfloat_remStepMBy32(uint8_t size_words,
                       uint32_t const* remPtr,
                       uint8_t dist,
                       uint32_t const* bPtr,
                       uint32_t q,
                       uint32_t* zPtr)
{
    auto index = indexWordLo(size_words);
    auto const lastIndex = indexWordHi(size_words);
    uint64_t dwordProd = static_cast<uint64_t>(bPtr[index]) * q;
    uint32_t wordRem = remPtr[index];
    uint32_t wordShiftedRem = wordRem << dist;
    uint32_t wordProd = static_cast<uint32_t>(dwordProd);
    zPtr[index] = wordShiftedRem - wordProd;

    if (index != lastIndex) {
        uint8_t const uNegDist = 31u & -dist;
        bool borrow = wordShiftedRem < wordProd;

        for (;;) {
            wordShiftedRem = wordRem >> uNegDist;
            index += wordIncr;
            dwordProd = static_cast<uint64_t>(bPtr[index]) * q + (dwordProd >> 32);
            wordRem = remPtr[index];
            wordShiftedRem |= wordRem << dist;
            wordProd = static_cast<uint32_t>(dwordProd);
            zPtr[index] = wordShiftedRem - wordProd - !!(borrow);

            if (index == lastIndex) {
                break;
            }

            borrow = borrow ? wordShiftedRem <= wordProd : wordShiftedRem < wordProd;
        }
    }

}

}  // namespace internals
}  // namespace softfloat
