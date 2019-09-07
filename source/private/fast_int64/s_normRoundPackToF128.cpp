
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

namespace softfloat {
namespace internals {
namespace fast_int64 {

float128_t
softfloat_normRoundPackToF128(bool const sign,
                              int32_t exp,
                              uint64_t sig64,
                              uint64_t sig0)
{
    if (0 == sig64) {
        exp -= 64;
        sig64 = sig0;
        sig0 = 0;
    }

    auto const shiftDist = count_leading_zeros(sig64) - 15;
    exp -= shiftDist;

    if (0 <= shiftDist) {
        if (0 != shiftDist) {
            uint128 const sig128 = softfloat_shortShiftLeft128(uint128{sig64, sig0}, static_cast<uint8_t>(shiftDist));
            sig64 = sig128.v64;
            sig0 = sig128.v0;
        }

        if (static_cast<uint32_t>(exp) < 0x7FFD) {
            return static_cast<float128_t>(uint128{packToF128UI64(sign, 0 != (sig64 | sig0) ? exp : 0, sig64), sig0});
        }

        return softfloat_roundPackToF128(sign, exp, sig64, sig0, 0);
    }

    uint128_extra const sig128Extra = softfloat_shortShiftRightJam128Extra(uint128{sig64, sig0}, 0, static_cast<uint8_t>(-shiftDist));
    return softfloat_roundPackToF128(sign, exp, sig128Extra.v.v64, sig128Extra.v.v0, sig128Extra.extra);
}

}  // namespace fast_int64
}  // namespace internals
}  // namespace softfloat
