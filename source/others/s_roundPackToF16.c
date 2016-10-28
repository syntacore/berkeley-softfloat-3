
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

#include "internals.h"

#include "softfloat/functions.h"

float16_t
softfloat_roundPackToF16(bool sign, int16_t exp, uint16_t sig)
{
    uint8_t roundBits;

    uint8_t const roundingMode = softfloat_roundingMode;
    bool const roundNearEven = roundingMode == softfloat_round_near_even;
    uint8_t roundIncrement = 0x8;
    if (!roundNearEven && (roundingMode != softfloat_round_near_maxMag)) {
        roundIncrement =
            roundingMode == (sign ? softfloat_round_min : softfloat_round_max) ? 0xF : 0;
    }
    roundBits = sig & 0xF;
    if (0x1D <= (unsigned int)exp) {
        if (exp < 0) {
            bool const isTiny =
                softfloat_detectTininess == softfloat_tininess_beforeRounding ||
                exp < -1 ||
                sig + roundIncrement < 0x8000;
            sig = softfloat_shiftRightJam32(sig, -exp);
            exp = 0;
            roundBits = sig & 0xF;
            if (isTiny && roundBits) {
                softfloat_raiseFlags(softfloat_flag_underflow);
            }
        } else if ((0x1D < exp) || (0x8000 <= sig + roundIncrement)) {
            softfloat_raiseFlags(softfloat_flag_overflow | softfloat_flag_inexact);
            return ui16_as_f16(packToF16UI(sign, 0x1F, 0) - !roundIncrement);
        }
    }
    if (roundBits) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }
    sig = (sig + roundIncrement) >> 4;
    sig &= ~(uint16_t)(!(roundBits ^ 8) & roundNearEven);
    return ui16_as_f16(packToF16UI(sign, sig ? exp : 0, sig));
}
