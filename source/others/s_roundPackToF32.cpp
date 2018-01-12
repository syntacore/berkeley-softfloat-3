
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

#include "internals.hpp"
#include "softfloat/functions.h"
#include <assert.h>

/**
@param[in] sign sign bit
@param[in] exp exponent 8-bit field
@param[in] sig mantissa bits [30..0]
*/
float32_t
softfloat_roundPackToF32(bool sign, int16_t exp, uint32_t sig)
{
    assert(softfloat_round_near_even <= softfloat_roundingMode && softfloat_roundingMode <= softfloat_round_near_maxMag);
    enum softfloat_round_mode const roundingMode = (enum softfloat_round_mode)softfloat_roundingMode;
    uint8_t const roundIncrement =
        softfloat_round_near_even == roundingMode || softfloat_round_near_maxMag == roundingMode ? /* one half */ 0x40 :
        (sign ? softfloat_round_min : softfloat_round_max) == roundingMode ? 0x7F : 0;
    uint32_t roundBits = sig & ~(~UINT32_C(0) << 7);
    if (0xFD <= (uint16_t)exp) {
        if (exp < 0) {
            bool const isTiny =
                softfloat_detectTininess == softfloat_tininess_beforeRounding ||
                exp < -1 ||
                sig + roundIncrement < 0x80000000;
            sig = softfloat_shiftRightJam32(sig, -exp);
            exp = 0;
            roundBits = sig & ~(~UINT32_C(0) << 7);
            if (isTiny && 0 != roundBits) {
                softfloat_raiseFlags(softfloat_flag_underflow);
            }
        } else if (0xFD < exp || 0x80000000 <= sig + roundIncrement) {
            softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);
            return u_as_f_32(packToF32UI(sign, 0xFF, 0) - !roundIncrement);
        }
    }
    if (0 != roundBits) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }
    bool const exact_mid_point = 0 == (roundBits ^ 0x40);
    bool const to_even = exact_mid_point && softfloat_round_near_even == roundingMode;
    uint32_t const sig1 = ((sig + roundIncrement) >> 7) & ~(uint32_t)to_even;
    return u_as_f_32(packToF32UI(sign, 0 != sig1 ? exp : 0, sig1));
}
