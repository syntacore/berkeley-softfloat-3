
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

#include "softfloat/functions.h"

#include "internals.hpp"
#include "specialize.hpp"

uint32_t
f128_to_ui32_r_minMag(float128_t a,
                      bool exact)
{
    using namespace softfloat::internals;
    uint64_t const uiA64 = f_as_u_128(a).v64;
    uint64_t const uiA0 = f_as_u_128(a).v0;
    int32_t const exp = expF128UI64(uiA64);
    uint64_t sig64 = fracF128UI64(uiA64) | (uiA0 != 0);
    int32_t const shiftDist = 0x402F - exp;

    if (49 <= shiftDist) {
        if (exact && (exp | sig64)) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return 0;
    }

    bool const sign = signF128UI64(uiA64);

    if (sign || shiftDist < 17) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            exp == 0x7FFF && sig64 ? ui32_fromNaN :
            sign ? ui32_fromNegOverflow : ui32_fromPosOverflow;
    }

    sig64 |= UINT64_C(0x0001000000000000);
    uint32_t const z = static_cast<uint32_t>(sig64 >> shiftDist);

    if (exact && static_cast<uint64_t>(z) << shiftDist != sig64) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return z;
}

