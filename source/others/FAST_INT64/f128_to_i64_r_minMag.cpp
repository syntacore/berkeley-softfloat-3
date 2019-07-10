
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

int64_t
f128_to_i64_r_minMag(float128_t a,
                     bool exact)
{
    using namespace softfloat::internals;
    uint128 const uA{a};
    uint64_t const uiA64 = uA.v64;
    uint64_t const uiA0 = uA.v0;
    bool const sign = is_sign(uiA64);
    int32_t const exp = expF128UI64(uiA64);
    uint64_t sig64 = fracF128UI64(uiA64);
    uint64_t const sig0 = uiA0;

    int32_t const shiftDist = 0x402F - exp;

    if (shiftDist < 0) {
        if (shiftDist < -14) {
            if (UINT64_C(0xC03E000000000000) == uiA64 && (sig0 < UINT64_C(0x0002000000000000))) {
                if (exact && sig0) {
                    softfloat_raiseFlags(softfloat_flag_inexact);
                }

                return INT64_MIN;
            } else {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return
                    0x7FFF == exp && 0 != (sig64 | sig0) ? i64_fromNaN :
                    sign ? i64_fromNegOverflow : i64_fromPosOverflow;
            }
        } else {
            sig64 |= UINT64_C(0x0001000000000000);
            /** @todo Warning   C4244   '=': conversion from 'int32_t' to 'int8_t', possible loss of data */
            auto const negShiftDist = -shiftDist;
            int64_t const absZ = static_cast<int64_t>(sig64 << negShiftDist | sig0 >> (shiftDist & 63u));

            if (exact && static_cast<uint64_t>(sig0 << negShiftDist)) {
                softfloat_raiseFlags(softfloat_flag_inexact);
            }

            return sign ? -absZ : absZ;
        }
    } else {
        if (49 <= shiftDist) {
            if (exact && (exp | sig64 | sig0)) {
                softfloat_raiseFlags(softfloat_flag_inexact);
            }

            return 0;
        } else {
            sig64 |= UINT64_C(0x0001000000000000);
            int64_t const absZ = static_cast<int64_t>(sig64 >> shiftDist);

            if (exact && (sig0 || static_cast<uint64_t>(absZ << shiftDist) != sig64)) {
                softfloat_raiseFlags(softfloat_flag_inexact);
            }

            return sign ? -absZ : absZ;
        }
    }
}
