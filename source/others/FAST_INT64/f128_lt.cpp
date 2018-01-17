
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
All rights reserved.

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

bool
f128_lt(float128_t a,
        float128_t b)
{
    using namespace softfloat;
    uint64_t const uiA64 = reinterpret_cast<uint128 const&>(a).v64;
    uint64_t const uiA0 = reinterpret_cast<uint128 const&>(a).v0;
    uint64_t const uiB64 = reinterpret_cast<uint128 const&>(b).v64;
    uint64_t const uiB0 = reinterpret_cast<uint128 const&>(b).v0;

    if (isNaNF128UI(uiA64, uiA0) || isNaNF128UI(uiB64, uiB0)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return false;
    }

    bool const signA = signF128UI64(uiA64);
    bool const signB = signF128UI64(uiB64);
    return
        signA != signB ? signA && 0 != (((uiA64 | uiB64) & UINT64_C(0x7FFFFFFFFFFFFFFF)) | uiA0 | uiB0) :
        (uiA64 != uiB64 || uiA0 != uiB0) && signA != softfloat_lt128(uiA64, uiA0, uiB64, uiB0);
}

