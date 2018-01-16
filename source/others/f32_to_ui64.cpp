
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
#include "specialize.hpp"

/** @todo split to different implementations */
uint64_t
f32_to_ui64(float32_t a,
            uint8_t roundingMode,
            bool exact)
{
    using namespace softfloat;
    uint32_t const uiA = f_as_u_32(a);
    bool const sign = signF32UI(uiA);
    int16_t const exp = expF32UI(uiA);
    uint32_t sig = fracF32UI(uiA);
    int16_t const shiftDist = 0xBE - exp;

    if (shiftDist < 0) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            exp == 0xFF && sig ? ui64_fromNaN :
            sign ? ui64_fromNegOverflow : ui64_fromPosOverflow;
    }

    if (exp) {
        sig |= 0x00800000;
    }

#ifdef SOFTFLOAT_FAST_INT64
    uint64_t sig64 = (uint64_t)sig << 40;
    uint64_t extra = 0;

    if (shiftDist) {
        uint64_extra const sig64Extra = softfloat_shiftRightJam64Extra(sig64, 0, shiftDist);
        sig64 = sig64Extra.v;
        extra = sig64Extra.extra;
    }

    return
        softfloat_roundPackToUI64(sign, sig64, extra, roundingMode, exact);
#else
    uint32_t extSig[3];
    extSig[indexWord(3, 2)] = sig << 8;
    extSig[indexWord(3, 1)] = 0;
    extSig[indexWord(3, 0)] = 0;

    if (shiftDist) {
        softfloat_shiftRightJam96M(extSig, static_cast<uint8_t>(shiftDist), extSig);
    }

    return softfloat_roundPackMToUI64(sign, extSig, roundingMode, exact);
#endif
}

