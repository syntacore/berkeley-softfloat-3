
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

/**
@todo refactoring
*/
int64_t
f128_to_i64(float128_t const a,
            softfloat_round_mode const roundingMode,
            bool const exact)
{
    using namespace softfloat::internals::fast_int64;

    uint128 const uA{a};
    bool const sign = is_sign(uA.v64);
    int32_t const exp = exp_F128_UI64(uA.v64);
    uint64_t sig64 = frac_F128_UI64(uA.v64);

    int32_t const shiftDist = 0x402F - exp;

    if (shiftDist <= 0) {
        if (shiftDist < -15) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return
                0x7FFF == exp && 0 != (sig64 | uA.v0) ?
                i64_fromNaN :
                sign ?
                i64_fromNegOverflow :
                i64_fromPosOverflow;
        }

        sig64 |= UINT64_C(0x0001000000000000);

        if (shiftDist) {
            uint128 const sig128 = short_shift_left_128(uint128{sig64, uA.v0}, static_cast<uint8_t>(-shiftDist));
            return round_pack_to<int64_t>(sign, sig128.v64, sig128.v0, roundingMode, exact);
        }

        return round_pack_to<int64_t>(sign, sig64, uA.v0, roundingMode, exact);
    }

    if (0 != exp) {
        sig64 |= UINT64_C(0x0001000000000000);
    }

    uint64_extra const sigExtra = shift_right_jam_64Extra(sig64, uA.v0, static_cast<uint32_t>(shiftDist));
    return round_pack_to<int64_t>(sign, sigExtra.v, sigExtra.extra, roundingMode, exact);
}
