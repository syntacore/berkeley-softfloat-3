
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

int32_t
f64_to_i32_r_minMag(float64_t a,
                    bool exact)
{
    using namespace softfloat::internals;
    uint64_t const uiA = f_as_u_64(a);
    int16_t const exp = expF64UI(uiA);
    uint64_t const sig = fracF64UI(uiA);
    int16_t const shiftDist = 0x433 - exp;

    if (53 <= shiftDist) {
        if (exact && (exp | sig)) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return 0;
    }

    bool const sign = is_sign(uiA);

    if (shiftDist < 22) {
        if (sign && 0x41E == exp && sig < UINT64_C(0x0000000000200000)) {
            if (exact && 0 != sig) {
                softfloat_raiseFlags(softfloat_flag_inexact);
            }

            return INT32_MIN;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            0x7FF == exp && sig ? i32_fromNaN :
            sign ? i32_fromNegOverflow : i32_fromPosOverflow;
    }

    auto const sig_1 = sig | UINT64_C(0x0010000000000000);
    int32_t const absZ = static_cast<int32_t>(sig_1 >> shiftDist);

    if (exact && (static_cast<uint64_t>(static_cast<uint32_t>(absZ)) << shiftDist != sig_1)) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return sign ? -absZ : absZ;
}
