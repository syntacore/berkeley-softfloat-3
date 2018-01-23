
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
#include "target.hpp"

uint64_t
f32_to_ui64_r_minMag(float32_t a,
                     bool exact)
{
    using namespace softfloat::internals;
    uint32_t const uiA = f_as_u_32(a);
    int16_t const exp = expF32UI(uiA);
    uint32_t sig = fracF32UI(uiA);
    int16_t shiftDist = 0xBE - exp;

    if (64 <= shiftDist) {
        if (exact && (exp | sig)) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return 0;
    }

    bool const sign = signF32UI(uiA);

    if (sign || (shiftDist < 0)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            exp == 0xFF && sig ? ui64_fromNaN :
            sign ? ui64_fromNegOverflow : ui64_fromPosOverflow;
    }

    uint64_t sig64;
    uint64_t z;
    sig |= 0x00800000;
    sig64 = static_cast<uint64_t>(sig) << 40;
    z = sig64 >> shiftDist;
    shiftDist = 40 - shiftDist;

    if (exact && (shiftDist < 0) && static_cast<uint32_t>(sig << (shiftDist & 31))) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return z;
}

