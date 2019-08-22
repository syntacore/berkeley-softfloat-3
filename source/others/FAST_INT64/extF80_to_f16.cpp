
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

#ifndef SOFTFLOAT_FAST_INT64
#error Fast int64_t operations only
#endif

float16_t
extF80_to_f16(extFloat80_t a)
{
    using namespace softfloat::internals;
    uint16_t const uiA64 = a.signExp;
    uint64_t const uiA0 = a.signif;
    bool const sign = is_sign(uiA64);
    int32_t exp = expExtF80UI64(uiA64);
    uint64_t const sig = uiA0;

    if (exp == 0x7FFF) {
        if (sig & INT64_MAX) {
            return u_as_f(softfloat_commonNaNToF16UI(softfloat_extF80UIToCommonNaN(uiA64, uiA0)));
        }

        return u_as_f(packToF16UI(sign, 0x1F, 0));
    }

    uint16_t const sig16 = static_cast<uint16_t>(softfloat_shortShiftRightJam64(sig, 49));

    if (!(exp | sig16)) {
        return u_as_f(packToF16UI(sign, 0, 0));
    }

    exp -= 0x3FF1;

    if (exp < -0x40) {
        exp = -0x40;
    }

    /**
    @todo Warning   C4242   'function': conversion from 'int32_t' to 'int16_t', possible loss of data
    */
    return softfloat_roundPackToF16(sign, static_cast<int16_t>(exp), sig16);
}
