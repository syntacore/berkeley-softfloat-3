
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All rights reserved.
*/
/*
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
f32_sqrt(float32_t a)
{
    using namespace softfloat::internals;
    uint32_t const uiA = f_as_u(a);
    bool const signA = is_sign(uiA);
    int16_t expA = get_exp(uiA);
    uint32_t sigA = get_frac(uiA);

    if (0xFF == expA) {
        if (0 != sigA) {
            return u_as_f_32(propagate_NaN(uiA, UINT32_C(0)));
        }

        if (!signA) {
            return a;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_32(defaultNaNF32UI);
    }

    if (signA) {
        if (0 == (expA | sigA)) {
            return a;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_32(defaultNaNF32UI);
    }

    if (0 == expA) {
        if (0 == sigA) {
            return a;
        }

        exp16_sig32 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int16_t const expZ = ((expA - 0x7F) >> 1) + 0x7E;
    expA &= 1;
    sigA = (sigA | 0x00800000) << 8;
    uint32_t sigZ = (static_cast<uint64_t>(sigA) * softfloat_approxRecipSqrt32_1(static_cast<uint32_t>(expA), sigA)) >> 32;

    if (expA) {
        sigZ >>= 1;
    }

    sigZ += 2;

    if ((sigZ & 0x3F) < 2) {
        uint32_t const shiftedSigZ = sigZ >> 2;
        uint32_t const negRem = shiftedSigZ * shiftedSigZ;
        sigZ &= ~3;

        if (negRem & 0x80000000) {
            sigZ |= 1;
        } else if (negRem) {
            --sigZ;
        }
    }

    return softfloat_roundPackToF32(0, expZ, sigZ);
}
