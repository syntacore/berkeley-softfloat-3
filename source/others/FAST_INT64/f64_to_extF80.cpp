
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

#include "target.hpp"

extFloat80_t
f64_to_extF80(float64_t a)
{
    using namespace softfloat::internals;
    uint128 uiZ;
    uint16_t uiZ64;
    uint64_t uiZ0;


    uint64_t const uiA = f_as_u_64(a);
    bool const sign = signF64UI(uiA);
    int16_t exp = expF64UI(uiA);
    uint64_t frac = fracF64UI(uiA);

    if (exp == 0x7FF) {
        if (frac) {
            uiZ = softfloat_commonNaNToExtF80UI(softfloat_f64UIToCommonNaN(uiA));
            uiZ64 = static_cast<uint16_t>(uiZ.v64);
            uiZ0 = uiZ.v0;
        }
        else {
            uiZ64 = packToExtF80UI64(sign, 0x7FFF);
            uiZ0 = UINT64_C(0x8000000000000000);
        }

        goto uiZ;
    }

    if (!exp) {
        if (!frac) {
            uiZ64 = packToExtF80UI64(sign, 0);
            uiZ0 = 0;
            goto uiZ;
        }

        exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(frac);
        exp = normExpSig.exp;
        frac = normExpSig.sig;
    }

    uiZ64 = packToExtF80UI64(sign, exp + 0x3C00u);
    uiZ0 = (frac | UINT64_C(0x0010000000000000)) << 11;
uiZ:
    extFloat80_t uZ;
    uZ.signExp = uiZ64;
    uZ.signif = uiZ0;
    return uZ;
}
