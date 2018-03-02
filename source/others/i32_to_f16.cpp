
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

float16_t
i32_to_f16(int32_t a)
{
    using namespace softfloat::internals;
    bool const sign = a < 0;
    /** @bug INT32_MIN */
    uint32_t absA = static_cast<uint32_t>(sign ? -a : a);
    int8_t const shiftDist = softfloat_countLeadingZeros32(absA) - 21;

    if (0 <= shiftDist) {
        return u_as_f_16(a ? packToF16UI(sign, 0x18 - shiftDist, static_cast<uint16_t>(absA << shiftDist)) : 0u);
    }

    int8_t const shiftDist1 = shiftDist + 4;
    uint16_t const sig =
        static_cast<uint16_t>(
            shiftDist1 < 0 ? absA >> -shiftDist1 | !!(0 != (absA << (shiftDist1 & 0x1F))) :
            absA << shiftDist1);
    return softfloat_roundPackToF16(sign, 0x1C - shiftDist1, sig);
}
