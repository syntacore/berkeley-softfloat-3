
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
f32_to_i64_r_minMag(float32_t const a,
                    bool const exact)
{
    using namespace softfloat::internals;
    int16_t const exp = get_exp(a);
    uint32_t const sig = get_frac(a);
    int16_t const shiftDist = 0xBE - exp;

    if (64 <= shiftDist) {
        if (exact && (exp | sig)) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return 0;
    }

    bool const sign = is_sign(a);

    if (shiftDist <= 0) {
        if (f_as_u(a) == pack_to_F32_UI(true, 0xBE, 0)) {
            return INT64_MIN;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            is_NaN(a) ? i64_fromNaN :
            sign ? i64_fromNegOverflow : i64_fromPosOverflow;
    }

    auto const sig_1 = sig | 0x00800000;
    uint64_t const sig64 = static_cast<uint64_t>(sig_1) << 40;
    int64_t const absZ = static_cast<int64_t>(sig64 >> shiftDist);
    auto const shiftDist1 = 40 - shiftDist;

    if (exact && shiftDist1 < 0 && 0 != static_cast<uint32_t>(sig_1 << (shiftDist1 & 31))) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return sign ? -absZ : absZ;
}
