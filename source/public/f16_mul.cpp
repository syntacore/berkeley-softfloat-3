
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

float16_t
f16_mul(float16_t const a,
        float16_t const b)
{
    using namespace softfloat::internals;

    bool const signA = is_sign(a);
    bool const signB = is_sign(b);
    bool const signZ = signA != signB;

    if (is_NaN(a) || is_NaN(b)) {
        return propagate_NaN(a, b);
    }

    if (is_inf(a) || is_inf(b)) {
        if (is_zero(a) || is_zero(b)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f(defaultNaNF16UI);
        }

        return make_signed_inf<float16_t>(signZ);
    }

    if (is_zero(a) || is_zero(b)) {
        return make_signed_zero<float16_t>(signZ);
    }

    int8_t expA = get_exp(a);
    uint16_t sigA = get_frac(a);

    if (0 == expA) {
        exp8_sig16 const normExpSig{sigA};
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int8_t expB = get_exp(b);
    uint16_t sigB = get_frac(b);

    if (0 == expB) {
        exp8_sig16 const normExpSig{sigB};
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int8_t expZ = expA + expB - 0xF;
    sigA = static_cast<uint16_t>((sigA | 0x0400) << 4);
    sigB = static_cast<uint16_t>((sigB | 0x0400) << 5);
    uint32_t const sig32Z = static_cast<uint32_t>(sigA) * sigB;
    uint16_t sigZ = sig32Z >> 16;

    if (0 != (sig32Z & 0xFFFF)) {
        sigZ |= 1;
    }

    if (sigZ < 0x4000) {
        --expZ;
        sigZ <<= 1;
    }

    return round_pack_to_F16(signZ, expZ, sigZ);
}
