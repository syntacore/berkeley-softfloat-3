
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

#include "model.hpp"

void
f128M_to_extF80M(const float128_t* aPtr,
                 extFloat80_t* zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    *zPtr = f128_to_extF80(*aPtr);
#else
    using namespace softfloat::internals::slow_int64;
    const uint32_t* aWPtr;
    extFloat80M* zSPtr;
    uint32_t uiA96;
    bool sign;
    int32_t exp;
    uint32_t sig[4];


    aWPtr = (const uint32_t*)aPtr;
    zSPtr = zPtr;

    uiA96 = aWPtr[index_word_hi(4)];
    sign = is_sign(uiA96);
    exp = exp_F128_UI96(uiA96);

    if (exp == 0x7FFF) {
        if (is_NaN_M_F128(aWPtr)) {
            *zSPtr = extFloat80M(commonNaN_from_M_f128(aWPtr));
            return;
        }

        zSPtr->signExp = pack_to_extF80_UI64(sign, 0x7FFF);
        zSPtr->signif = UINT64_C(0x8000000000000000);
        return;
    }

    exp = shift_norm_sig_M_F128(aWPtr, 15, sig);

    if (exp == -128) {
        zSPtr->signExp = pack_to_extF80_UI64(sign, 0);
        zSPtr->signif = 0;
        return;
    }

    if (sig[index_word(4, 0)]) {
        sig[index_word(4, 1)] |= 1;
    }

    round_pack_to_M_extF80M(
        sign, exp, &sig[index_multiword_hi(4, 3)], 80, zSPtr);

#endif
}
