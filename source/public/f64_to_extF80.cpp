
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

extFloat80_t
f64_to_extF80(float64_t const a)
{
    using namespace softfloat::internals;
    bool const sign = is_sign(a);
    int16_t exp = get_exp(a);
    uint64_t frac = get_frac(a);

    if (!is_finite(a)) {
        if (is_NaN(a)) {
            uint128 const uiZ = uint128(commonNaN_from_f64UI(f_as_u(a)));
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(uiZ.v64);
            uZ.signif = uiZ.v0;
            return uZ;
        } else {
            extFloat80_t uZ;
            uZ.signExp = pack_to_extF80_UI64(sign, 0x7FFF);
            uZ.signif = UINT64_C(0x8000000000000000);
            return uZ;
        }
    }

    if (is_zero(a)) {
        extFloat80_t uZ;
        uZ.signExp = pack_to_extF80_UI64(sign, 0);
        uZ.signif = 0;
        return uZ;
    }

    if (0 == exp) {
        exp16_sig64 const normExpSig(frac);
        exp = normExpSig.exp;
        frac = normExpSig.sig;
    }

    extFloat80_t uZ;
    uZ.signExp = pack_to_extF80_UI64(sign, exp + 0x3C00u);
    uZ.signif = (frac | UINT64_C(0x0010000000000000)) << 11;
    return uZ;
}
