
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All rights reserved.
*/
/*
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

uint32_t
f16_to_ui32(float16_t const a,
            softfloat_round_mode const roundingMode,
            bool const exact)
{
    using namespace softfloat::internals;

    bool const sign = is_sign(a);
    uint16_t const frac = get_frac(a);

    if (!is_finite(a)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            frac ? ui32_fromNaN :
            sign ? ui32_fromNegOverflow :
            ui32_fromPosOverflow;
    }

    uint32_t sig32 = frac;
    int8_t const exp = get_exp(a);

    if (0 != exp) {
        sig32 |= 0x0400;
        {
            int8_t const shiftDist = exp - 0x19;

            if (0 <= shiftDist && !sign) {
                return sig32 << shiftDist;
            }
        }
        int8_t const shiftDist = exp - 0x0D;

        if (0 < shiftDist) {
            sig32 <<= shiftDist;
        }
    }

    return round_pack_to<uint32_t>(sign, sig32, roundingMode, exact);
}
