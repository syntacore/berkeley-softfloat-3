
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

#include "internals.hpp"
#include "specialize.hpp"
#include "softfloat/functions.h"

namespace softfloat {
namespace internals {

int64_t
softfloat_roundPackMToI64(bool sign,
                          uint32_t* extSigPtr,
                          uint8_t roundingMode,
                          bool exact)
{
    int64_t z;


    bool const roundNearEven = (roundingMode == softfloat_round_near_even);
    uint32_t const sigExtra = extSigPtr[indexWordLo(3)];
    bool doIncrement =
        !roundNearEven &&  softfloat_round_near_maxMag != roundingMode ? (roundingMode == (sign ? softfloat_round_min : softfloat_round_max)) && sigExtra :
        0x80000000 <= sigExtra;
    uint64_t sig = static_cast<uint64_t>(extSigPtr[indexWord(3, 2)]) << 32 | extSigPtr[indexWord(3, 1)];

    if (doIncrement) {
        ++sig;

        if (!sig) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return sign ? i64_fromNegOverflow : i64_fromPosOverflow;
        }

        if (!(sigExtra & 0x7FFFFFFF) && roundNearEven) {
            sig &= ~1;
        }
    }

    z = sign ? -static_cast<int64_t>(sig) : static_cast<int64_t>(sig);

    if (z && ((z < 0) ^ sign)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return sign ? i64_fromNegOverflow : i64_fromPosOverflow;
    }

    if (exact && sigExtra) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return z;
}

}  // namespace internals
}  // namespace softfloat
