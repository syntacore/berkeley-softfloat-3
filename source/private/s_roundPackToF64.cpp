
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015 The Regents of the University of
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

namespace softfloat {
namespace internals {

float64_t
round_pack_to_F64(bool sign,
                  int16_t exp,
                  uint64_t sig)
{
    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
    uint16_t const roundIncrement =
        softfloat_round_near_even == softfloat_roundingMode || softfloat_round_near_maxMag == softfloat_roundingMode ? 0x200u :
        (sign ? softfloat_round_min : softfloat_round_max) == softfloat_roundingMode ? 0x3FFu :
        0u;
    uint16_t roundBits = sig & 0x3FF;

    if (0x7FDu <= static_cast<uint16_t>(exp)) {
        if (exp < 0) {
            bool const isTiny =
                softfloat_tininess_beforeRounding == detect_tininess ||
                exp < -1 ||
                sig + roundIncrement < uint64_t(INT64_MIN);
            sig = shift_right_jam_64(sig, static_cast<uint32_t>(-exp));
            exp = 0;
            roundBits = sig & 0x3FF;

            if (isTiny && roundBits) {
                softfloat_raiseFlags(softfloat_flag_underflow);
            }
        } else if (0x7FD < exp || uint64_t(INT64_MIN) <= sig + roundIncrement) {
            softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);
            return u_as_f(pack_to_F64_UI(sign, 0x7FF, 0) - !roundIncrement);
        }
    }

    if (0 != roundBits) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    sig = ((sig + roundIncrement) >> 10) & (~static_cast<uint64_t>(0 == (roundBits ^ 0x200) && softfloat_round_near_even == softfloat_roundingMode));
    return u_as_f(pack_to_F64_UI(sign, 0 != sig ? exp : 0, sig));
}

}  // namespace internals
}  // namespace softfloat
