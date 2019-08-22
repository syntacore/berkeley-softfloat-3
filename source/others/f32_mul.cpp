
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

float32_t
f32_mul(float32_t a,
        float32_t b)
{
    using namespace softfloat::internals;
    uint32_t const uiA = f_as_u(a);
    uint32_t const uiB = f_as_u(b);

    if (is_NaN(uiA) || is_NaN(uiB)) {
        return u_as_f(propagate_NaN(uiA, uiB));
    }

    bool const signA = is_sign(uiA);
    int16_t expA = get_exp(uiA);
    uint32_t sigA = get_frac(uiA);

    bool const signB = is_sign(uiB);
    int16_t expB = get_exp(uiB);
    uint32_t sigB = get_frac(uiB);

    bool const signZ = signA != signB;

    if (expA == 0xFF || expB == 0xFF) {
        bool const is_undefined =
            /* a is infinity */ expA == 0xFF ? /* b is zero */ 0 == expB && 0 == sigB :
            /* b is infinity and a is zero */ 0 == expA && 0 == sigA;

        if (is_undefined) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f(defaultNaNF32UI);
        }

        return make_signed_inf<float32_t>(signZ);
    }

    if (0 == expA) {
        if (0 == sigA) {
            return make_signed_zero<float32_t>(signZ);
        }

        exp16_sig32 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (0 == expB) {
        if (0 == sigB) {
            return make_signed_zero<float32_t>(signZ);
        }

        exp16_sig32 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int16_t expZ = expA + expB - 127;
    sigA = (sigA | 0x00800000) << 7;
    sigB = (sigB | 0x00800000) << 8;
    uint32_t sigZ = static_cast<uint32_t>(softfloat_shortShiftRightJam64(static_cast<uint64_t>(sigA) * sigB, 32));

    if (sigZ < 0x40000000) {
        --expZ;
        sigZ <<= 1;
    }

    return softfloat_roundPackToF32(signZ, expZ, sigZ);
}
