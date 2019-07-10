
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

uint32_t
f32_to_ui32(float32_t const a,
            softfloat_round_mode const roundingMode,
            bool const exact)
{
    using namespace softfloat::internals;
    uint32_t const uiA = f_as_u_32(a);
    bool const sign = signF32UI(uiA);
    int16_t const exp = expF32UI(uiA);

    if (isNaNF32UI(uiA)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return ui32_fromNaN;
    } else if (isZero32UI(uiA)) {
        return 0u;
    } else if (isInf32UI(uiA)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return sign ? ui32_fromNegOverflow : ui32_fromPosOverflow;
    } else {
        uint64_t const sig64 = static_cast<uint64_t>(fracF32UI(uiA) | (0 != exp? 0x00800000: 0) ) << 32;
        int16_t const shiftDist = 0xAA - exp;
        return
            roundPackTo<uint32_t>(sign,
                                  0 < shiftDist ?
                                  softfloat_shiftRightJam64(sig64, static_cast<uint32_t>(shiftDist)) :
                                  sig64,
                                  roundingMode,
                                  exact);
    }
}
