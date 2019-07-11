
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

int64_t
f64_to_i64_r_minMag(float64_t const a,
                    bool const exact)
{
    using namespace softfloat::internals;
    uint64_t const uiA = f_as_u_64(a);
    bool const sign = is_sign(uiA);
    int16_t const exp = expF64UI(uiA);
    uint64_t const sig = fracF64UI(uiA);
    int16_t const shiftDist = 0x433 - exp;

    if (shiftDist <= 0) {
        if (shiftDist < -10) {
            if (uiA == packToF64UI(1, 0x43E, 0)) {
                return INT64_MIN;
            }

            softfloat_raiseFlags(softfloat_flag_invalid);
            return
                exp == 0x7FF && sig ? i64_fromNaN :
                sign ? i64_fromNegOverflow : i64_fromPosOverflow;
        }

        int64_t const absZ = static_cast<int64_t>((sig | UINT64_C(0x0010000000000000)) << -shiftDist);
        return sign ? -absZ : absZ;
    }

    if (53 <= shiftDist) {
        if (exact && 0 != (exp | sig)) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return 0;
    }

    uint64_t sig1 = sig | UINT64_C(0x0010000000000000);
    uint64_t const absZ = sig1 >> shiftDist;

    if (exact && (absZ << shiftDist) != sig1) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return sign ? -static_cast<int64_t>(absZ) : static_cast<int64_t>(absZ);
}
