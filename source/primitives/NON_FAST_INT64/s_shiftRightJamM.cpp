
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

#include "primitives/functions.hpp"

#include <assert.h>
#include <string.h>

void
softfloat_shiftRightJamM(
    uint8_t size_words,
    const uint32_t *aPtr,
    uint32_t dist,
    uint32_t *zPtr
)
{
    uint32_t *ptr = NULL;

    uint32_t wordJam = 0;
    uint32_t wordDist = dist >> 5;
    if (wordDist) {
        if (size_words < wordDist) {
            wordDist = size_words;
        }
        ptr = (uint32_t *)(aPtr + indexMultiwordLo(size_words, wordDist));
        {
            int i = wordDist;
            do {
                wordJam = *ptr++;
                if (wordJam) {
                    break;
                }
            } while (--i);
        }
        ptr = zPtr;
    }
    if (wordDist < size_words) {
        uint8_t innerDist;
        aPtr += indexMultiwordHiBut(size_words, wordDist);
        innerDist = dist & 31;
        if (innerDist) {
            /** @todo Warning	C4244	'=': conversion from 'uint32_t' to 'uint8_t', possible loss of data */
            softfloat_shortShiftRightJamM(
                size_words - wordDist,
                aPtr,
                innerDist,
                zPtr + indexMultiwordLoBut(size_words, wordDist)
            );
            if (!wordDist) goto wordJam;
        } else {
            aPtr += indexWordLo(size_words - wordDist);
            ptr = zPtr + indexWordLo(size_words);
            for (int i = size_words - wordDist; i; --i) {
                *ptr = *aPtr;
                aPtr += wordIncr;
                ptr += wordIncr;
            }
        }
        ptr = zPtr + indexMultiwordHi(size_words, wordDist);
    }
    assert(ptr);
    memset(ptr, 0, wordDist * (sizeof *ptr));
wordJam:
    if (wordJam) {
        zPtr[indexWordLo(size_words)] |= 1;
    }
}
