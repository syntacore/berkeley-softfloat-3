
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

#include "softfloat/functions.h"

#include "internals.h"
#include "specialize.h"

/** @todo split to different implementations */
#ifdef SOFTFLOAT_FAST_INT64

void
f16_to_extF80M(float16_t a, extFloat80_t *zPtr)
{
    *zPtr = f16_to_extF80(a);
}

#else

void f16_to_extF80M(float16_t a, extFloat80_t *zPtr)
{
    uint16_t uiZ64;
    uint32_t uiZ32;

    /** @bug cast to same type */
    struct extFloat80M *const zSPtr = (struct extFloat80M *) zPtr;
    uint16_t const uiA = f_as_u_16(a);
    bool const sign = signF16UI(uiA);
    int8_t exp = expF16UI(uiA);
    uint16_t frac = fracF16UI(uiA);

    if (exp == 0x1F) {
        if (frac) {
            struct commonNaN commonNaN;
            softfloat_f16UIToCommonNaN(uiA, &commonNaN);
            softfloat_commonNaNToExtF80M(&commonNaN, zSPtr);
            return;
        } else {
            uiZ64 = packToExtF80UI64(sign, 0x7FFF);
            uiZ32 = 0x80000000;
            goto uiZ;
        }
    }
    if (!exp) {
        if (!frac) {
            uiZ64 = packToExtF80UI64(sign, 0);
            uiZ32 = 0;
            goto uiZ;
        } else {
            struct exp8_sig16 const normExpSig = softfloat_normSubnormalF16Sig(frac);
            exp = normExpSig.exp;
            frac = normExpSig.sig;
        }
    }

    uiZ64 = packToExtF80UI64(sign, exp + 0x3FF0);
    uiZ32 = 0x80000000 | (uint32_t)frac << 21;
uiZ:
    zSPtr->signExp = uiZ64;
    zSPtr->signif = (uint64_t)uiZ32 << 32;
}

#endif
