
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
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
f16_sqrt(float16_t a)
{
    using namespace softfloat::internals;
    bool const signA = is_sign(a);
    int8_t expA = get_exp(a);
    uint16_t sigA = get_frac(a);

    if (is_NaN(a)) {
        return propagate_NaN(a, a);
    }

    if (is_inf(a)) {
        if (!signA) {
            return a;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f(defaultNaNF16UI);
    }

    if (is_zero(a)) {
        return a;
    }

    if (signA) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f(defaultNaNF16UI);
    }

    if (0 == expA) {
        exp8_sig16 const normExpSig{sigA};
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int8_t expZ = ((expA - 0xF) >> 1) + 0xE;
    expA &= 1;
    sigA |= 0x0400;
    int const index = (sigA >> 6 & 0xE) + expA;
    uint16_t const r0 =
        softfloat_approxRecipSqrt_1k0s[index] -
        ((static_cast<uint32_t>(softfloat_approxRecipSqrt_1k1s[index]) * (sigA & 0x7F)) >> 11);
    uint32_t ESqrR0 = (static_cast<uint32_t>(r0) * r0) >> 1;

    if (expA) {
        ESqrR0 >>= 1;
    }

    uint16_t const sigma0 = static_cast<uint16_t>(~((ESqrR0 * sigA) >> 16));
    uint16_t recipSqrt16 = r0 + ((static_cast<uint32_t>(r0) * sigma0) >> 25);

    if (!(recipSqrt16 & 0x8000)) {
        recipSqrt16 = 0x8000;
    }

    uint16_t sigZ = (static_cast<uint32_t>(sigA << 5) * recipSqrt16) >> 16;

    if (expA) {
        sigZ >>= 1;
    }

    ++sigZ;

    if (0 == (sigZ & 7)) {
        uint16_t const shiftedSigZ = static_cast<uint16_t>(sigZ >> 1);
        uint16_t const negRem = static_cast<uint16_t>(shiftedSigZ * shiftedSigZ);
        sigZ &= ~1;

        if (negRem & 0x8000) {
            sigZ |= 1;
        } else {
            if (negRem) {
                --sigZ;
            }
        }
    }

    return softfloat_roundPackToF16(0, expZ, sigZ);
}
