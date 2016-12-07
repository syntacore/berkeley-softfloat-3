
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

#include "primitives/functions.h"

#include <assert.h>
#include <string.h>

static
void
softfloat_shortShiftRightJamM(
    uint8_t size_words,
    const uint64_t *aPtr,
    uint8_t dist,
    uint64_t *zPtr
)
{
    uint8_t const uNegDist = -(int8_t)dist;
    unsigned index = indexWordLo(size_words);
    unsigned const lastIndex = indexWordHi(size_words);
    uint64_t wordA = aPtr[index];
    uint64_t partWordZ = wordA >> dist;
    if (partWordZ << dist != wordA) {
        partWordZ |= 1;
    }
    while (index != lastIndex) {
        wordA = aPtr[index + wordIncr];
        zPtr[index] = wordA << (uNegDist & 63) | partWordZ;
        index += wordIncr;
        partWordZ = wordA >> dist;
    }
    zPtr[index] = partWordZ;

}

void
softfloat_shiftRightJam256M(
    const uint64_t *aPtr, uint32_t dist, uint64_t *zPtr)
{
    uint64_t *ptr = NULL;
    uint8_t innerDist;

    uint64_t wordJam = 0;
    uint32_t wordDist = dist >> 6;
    if (wordDist) {
        if (4 < wordDist) {
            wordDist = 4;
        }
        ptr = (uint64_t *)(aPtr + indexMultiwordLo(4, wordDist));
        {
            uint32_t i = wordDist;
            do {
                wordJam = *ptr++;
                if (wordJam) {
                    break;
                }
            } while (--i);
        }
        ptr = zPtr;
    }
    if (wordDist < 4) {
        aPtr += indexMultiwordHiBut(4, wordDist);
        innerDist = dist & 63;
        if (innerDist) {
            /** @todo Warning	C4244	'=': conversion from 'uint32_t' to 'uint8_t', possible loss of data */
            softfloat_shortShiftRightJamM(
                4 - wordDist,
                aPtr,
                innerDist,
                zPtr + indexMultiwordLoBut(4, wordDist)
            );
            if (!wordDist) {
                goto wordJam;
            }
        } else {
            aPtr += indexWordLo(4 - wordDist);
            ptr = zPtr + indexWordLo(4);
            for (int i = 4 - wordDist; i; --i) {
                *ptr = *aPtr;
                aPtr += wordIncr;
                ptr += wordIncr;
            }
        }
        ptr = zPtr + indexMultiwordHi(4, wordDist);
    }
    assert(ptr);
    memset(ptr, 0, wordDist * (sizeof *ptr));
wordJam:
    if (wordJam) {
        zPtr[indexWordLo(4)] |= 1;
    }
}
