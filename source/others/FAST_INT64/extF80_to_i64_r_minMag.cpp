
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

#include "softfloat/functions.h"

#include "internals.hpp"
#include "specialize.hpp"

int64_t
extF80_to_i64_r_minMag(extFloat80_t const a,
                       bool exact)
{
    using namespace softfloat;
    uint16_t const uiA64 = a.signExp;
    int32_t const exp = expExtF80UI64(uiA64);
    uint64_t const sig = a.signif;
    int32_t const shiftDist = 0x403E - exp;

    if (64 <= shiftDist) {
        if (exact && (exp | sig)) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return 0;
    }

    bool const sign = signExtF80UI64(uiA64);

    if (shiftDist <= 0) {
        if (uiA64 == packToExtF80UI64(1, 0x403E) && sig == UINT64_C(0x8000000000000000)) {
            return -INT64_C(0x7FFFFFFFFFFFFFFF) - 1;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            exp == 0x7FFF && 0 != (sig & UINT64_C(0x7FFFFFFFFFFFFFFF)) ? i64_fromNaN : 
            sign ? i64_fromNegOverflow : i64_fromPosOverflow;
    }

    int64_t const absZ = static_cast<int64_t>(sig >> shiftDist);

    if (exact && static_cast<uint64_t>(sig << (-shiftDist & 63))) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return sign ? -absZ : absZ;

}

