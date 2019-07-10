
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

#include "target.hpp"

bool
extF80_lt(extFloat80_t const a,
          extFloat80_t const b)
{
    using namespace softfloat::internals;
    auto const uiA64 = a.signExp;
    auto const uiA0  = a.signif;
    auto const uiB64 = b.signExp;
    auto const uiB0  = b.signif;

    if (isNaNExtF80UI(uiA64, uiA0) || isNaNExtF80UI(uiB64, uiB0)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return false;
    }

    bool const signA = signExtF80UI64(uiA64);
    bool const signB = signExtF80UI64(uiB64);
    typedef typename std::conditional<(sizeof uiA0 < sizeof uiA64), decltype(uiA64), decltype(uiA0)>::type largest_type;
    return
        signA != signB ?
        signA && 0 != static_cast<largest_type>(((uiA64 | uiB64) & 0x7FFF) | uiA0 | uiB0) :
        (uiA64 != uiB64 || uiA0 != uiB0) && softfloat_lt128(uiA64, uiA0, uiB64, uiB0) != signA;
}

