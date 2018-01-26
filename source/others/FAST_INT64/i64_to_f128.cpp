
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

float128_t
i64_to_f128(int64_t a)
{
    using namespace softfloat::internals;

    if (0 == a) {
        return static_cast<float128_t>(uint128{0, 0});
    } else {
        bool const sign = (a < 0);
        /** @bug for INT64_MIN */
        uint64_t const absA = static_cast<uint64_t>(sign ? -a : a);
        int8_t const shiftDist = softfloat_countLeadingZeros64(absA) + 49;

        if (64 <= shiftDist) {
            return static_cast<float128_t>(uint128{packToF128UI64(sign, 0x406E - shiftDist, absA << (shiftDist - 64)), 0});
        } else {
            uint128 const zSig = softfloat_shortShiftLeft128(0, absA, static_cast<uint8_t>(shiftDist));
            return static_cast<float128_t>(uint128{packToF128UI64(sign, 0x406E - shiftDist, zSig.v64), zSig.v0});
        }
    }
}

