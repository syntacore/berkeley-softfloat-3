
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

uint32_t
extF80_to_ui32(extFloat80_t const a,
               uint8_t const roundingMode,
               bool const exact)
{
    uint16_t const uiA64 = a.signExp;
    bool sign = signExtF80UI64(uiA64);
    int32_t const exp = expExtF80UI64(uiA64);
    uint64_t sig = a.signif;

    if (ui32_fromNaN != ui32_fromPosOverflow || ui32_fromNaN != ui32_fromNegOverflow) {
        if ((exp == 0x7FFF) && (sig & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
            if (ui32_fromNaN == ui32_fromPosOverflow) {
                sign = 0;
            } else if (ui32_fromNaN == ui32_fromNegOverflow) {
                sign = 1;
            } else {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return ui32_fromNaN;
            }
        }
    }

    int32_t shiftDist = 0x4032 - exp;

    if (shiftDist <= 0) {
        shiftDist = 1;
    }

    sig = softfloat_shiftRightJam64(sig, static_cast<uint32_t>(shiftDist));
    return softfloat_roundPackToUI32(sign, sig, roundingMode, exact);
}
