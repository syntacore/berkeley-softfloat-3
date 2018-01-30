
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

#include "target.hpp"

int32_t
f32_to_i32_r_minMag(float32_t const a,
                    bool const exact)
{
    using namespace softfloat::internals;
    uint32_t const uiA = f_as_u_32(a);
    int16_t const exp = expF32UI(uiA);
    uint32_t const sig = fracF32UI(uiA);
    int16_t const shiftDist = 0x9E - exp;

    if (32 <= shiftDist) {
        if (exact && (exp | sig)) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return 0;
    } else {
        bool const sign = signF32UI(uiA);

        if (shiftDist <= 0) {
            if (uiA == packToF32UI(true, 0x9E, 0)) {
                return INT32_MIN;
            } else {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return
                    0xFF == exp && sig ? i32_fromNaN :
                    sign ? i32_fromNegOverflow : i32_fromPosOverflow;
            }
        } else {
            auto const sig_1 = (sig | 0x00800000) << 8;
            int32_t const absZ = static_cast<int32_t>(sig_1 >> shiftDist);

            if (exact && (static_cast<uint32_t>(absZ) << shiftDist != sig_1)) {
                softfloat_raiseFlags(softfloat_flag_inexact);
            }

            return sign ? -absZ : absZ;
        }
    }
}
