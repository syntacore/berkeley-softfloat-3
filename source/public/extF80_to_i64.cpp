
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
extF80_to_i64(extFloat80_t const a,
              softfloat_round_mode const roundingMode,
              bool const exact)
{
    using namespace softfloat::internals::fast_int64;

    bool const sign = is_sign(a);
    int32_t const exp = get_exp(a);
    uint64_t sig = a.signif;
    int32_t const shiftDist = 0x403E - exp;

    uint64_t sigExtra;
    if (shiftDist <= 0) {

        if (shiftDist) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return
                (exp == INT16_MAX) && (sig & INT64_MAX) ? i64_fromNaN :
                sign ? i64_fromNegOverflow : i64_fromPosOverflow;
        }

        sigExtra = 0;
    } else {
        uint64_extra const sig64Extra = shift_right_jam_64Extra(sig, 0u, static_cast<uint32_t>(shiftDist));
        sig = sig64Extra.v;
        sigExtra = sig64Extra.extra;
    }

    return round_pack_to<int64_t>(sign, sig, sigExtra, roundingMode, exact);
}

