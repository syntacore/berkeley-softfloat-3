
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
f64_to_i32(float64_t const a,
           softfloat_round_mode const roundingMode,
           bool const exact)
{
    using namespace softfloat::internals;
    static bool const fromNaN_is_same_as_pos_overflow = i32_fromNaN == i32_fromPosOverflow;
    static bool const fromNaN_is_same_as_neg_overflow = i32_fromNaN == i32_fromNegOverflow;
    static bool const fromNaN_is_same_as_both_overflow = fromNaN_is_same_as_pos_overflow && fromNaN_is_same_as_neg_overflow;
    static bool const fromNaN_is_same_as_any_overflow = fromNaN_is_same_as_pos_overflow || fromNaN_is_same_as_neg_overflow;
    uint64_t const uiA = f_as_u_64(a);
    bool sign = is_sign(uiA);
    int16_t const exp = expF64UI(uiA);
    uint64_t sig = fracF64UI(uiA);

    // TODO: check and re-factor
    if (!fromNaN_is_same_as_both_overflow && exp == 0x7FF && 0 != sig) {
        if (!fromNaN_is_same_as_any_overflow) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return i32_fromNaN;
        }

        sign = fromNaN_is_same_as_neg_overflow;
    }

    if (0 != exp) {
        sig |= UINT64_C(0x0010000000000000);
    }

    int16_t const shiftDist = 0x427 - exp;
    return roundPackTo<int32_t>(sign, 0 < shiftDist ? softfloat_shiftRightJam64(sig, static_cast<uint32_t>(shiftDist)) : sig, roundingMode, exact);
}
