
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015 The Regents of the University of
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

void
f32_to_extF80M(float32_t a,
               extFloat80_t* zPtr)
{
#ifdef SOFTFLOAT_FAST_INT64
    *zPtr = f32_to_extF80(a);
#else
    using namespace softfloat::internals::slow_int64;

    extFloat80M* const zSPtr = zPtr;
    bool const sign = is_sign(a);
    int16_t exp = get_exp(a);
    uint32_t frac = get_frac(a);

    if (is_NaN(a)) {
        *zSPtr = extFloat80M(softfloat_f32UIToCommonNaN(f_as_u(a)));
        return;
    }

    if (is_inf(a)) {
        zSPtr->signExp = packToExtF80UI64(sign, 0x7FFF);
        zSPtr->signif = static_cast<uint64_t>(0x80000000) << 32;
        return;
    }

    if (is_zero(a)) {
        zSPtr->signExp = packToExtF80UI64(sign, 0);
        zSPtr->signif = 0;
        return;
    }

    if (0 == exp) {
        exp16_sig32 const normExpSig(frac);
        exp = normExpSig.exp;
        frac = normExpSig.sig;
    }

    zSPtr->signExp = packToExtF80UI64(sign, static_cast<uint16_t>(exp + 0x3F80));
    zSPtr->signif = static_cast<uint64_t>(0x80000000 | static_cast<uint32_t>(frac) << 8) << 32;
#endif
}
