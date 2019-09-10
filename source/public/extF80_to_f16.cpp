
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

float16_t
extF80_to_f16(extFloat80_t const a)
{
    using namespace softfloat::internals::fast_int64;

    bool const sign = is_sign(a);
    int32_t exp = get_exp(a);

    if (0x7FFF == exp) {
        if (0 != (a.signif & INT64_MAX)) {
            return u_as_f(commonNaN_to_F16UI(commonNaN_from_extF80UI(a.signExp, a.signif)));
        }

        return make_signed_inf<float16_t>(sign);
    }

    uint16_t const sig16 = static_cast<uint16_t>(short_shift_right_jam_64(a.signif, 49));

    if (0 == (exp | sig16)) {
        return u_as_f(pack_to_F16_UI(sign, 0, 0));
    }

    exp -= 0x3FF1;

    if (exp < -0x40) {
        exp = -0x40;
    }

    /**
    @todo Warning   C4242   'function': conversion from 'int32_t' to 'int16_t', possible loss of data
    */
    return round_pack_to_F16(sign, static_cast<int16_t>(exp), sig16);
}
