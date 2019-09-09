
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

int64_t
f32_to_i64(float32_t const a,
           softfloat_round_mode const roundingMode,
           bool const exact)
{
#if (SOFTFLOAT_FAST_INT64)
    using namespace softfloat::internals::fast_int64;

    bool const sign = is_sign(a);
    int16_t const exp = get_exp(a);
    uint32_t const sig = get_frac(a);

    auto const shiftDist = 0xBE - exp;

    if (shiftDist < 0) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            is_NaN(a) ? i64_fromNaN :
            sign ? i64_fromNegOverflow : i64_fromPosOverflow;
    }

    uint64_t const sig64 = static_cast<uint64_t>(sig | (0 != exp ? 0x00800000 : 0)) << 40;

    if (0 != shiftDist) {
        uint64_extra const sig64Extra = shift_right_jam_64Extra(sig64, 0u, static_cast<uint32_t>(shiftDist));
        return round_pack_to<int64_t>(sign, sig64Extra.v, sig64Extra.extra, roundingMode, exact);
    }

    return round_pack_to<int64_t>(sign, sig64, 0, roundingMode, exact);

#else
    using namespace softfloat::internals::slow_int64;

    bool const sign = is_sign(a);
    int16_t const exp = get_exp(a);
    uint32_t const sig = get_frac(a);

    int16_t const shiftDist = 0xBE - exp;

    if (shiftDist < 0) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            is_NaN(a) ? i64_fromNaN :
            sign ? i64_fromNegOverflow : i64_fromPosOverflow;
    }

    uint32_t extSig[3];
    extSig[index_word(3, 2)] = (sig | (0 != exp ? 0x00800000 : 0)) << 8;
    extSig[index_word(3, 1)] = 0;
    extSig[index_word(3, 0)] = 0;

    if (shiftDist) {
        shift_right_jam_M_96(extSig, static_cast<uint8_t>(shiftDist), extSig);
    }

    return round_pack_to_M<int64_t>(sign, extSig, roundingMode, exact);
#endif
}

