
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

#include "softfloat/functions.h"

#include "internals.hpp"
#include "specialize.hpp"

extern const uint16_t softfloat_approxRecipSqrt_1k0s[];
extern const uint16_t softfloat_approxRecipSqrt_1k1s[];

float16_t f16_sqrt(float16_t a)
{
    uint16_t const uiA = f_as_u_16(a);
    bool const signA = signF16UI(uiA);
    int8_t expA = expF16UI(uiA);
    uint16_t sigA = fracF16UI(uiA);

    if (expA == 0x1F) {
        if (sigA) {
            return u_as_f_16(softfloat_propagateNaNF16UI(uiA, 0));
        } else if (!signA) {
            return a;
        } else {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_16(defaultNaNF16UI);
        }
    } else if (signA) {
        if (!(expA | sigA)) {
            return a;
        } else {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_16(defaultNaNF16UI);
        }
    } else {
        if (!expA) {
            if (!sigA) {
                return a;
            } else {
                struct exp8_sig16 const normExpSig = softfloat_normSubnormalF16Sig(sigA);
                expA = normExpSig.exp;
                sigA = normExpSig.sig;
            }
        }

        int8_t expZ = ((expA - 0xF) >> 1) + 0xE;
        expA &= 1;
        sigA |= 0x0400;
        int const index = (sigA >> 6 & 0xE) + expA;
        uint16_t const r0 =
            softfloat_approxRecipSqrt_1k0s[index] - 
            (((uint32_t)softfloat_approxRecipSqrt_1k1s[index] * (sigA & 0x7F)) >> 11);
        uint32_t ESqrR0 = ((uint32_t)r0 * r0) >> 1;
        if (expA) {
            ESqrR0 >>= 1;
        }
        uint16_t const sigma0 = ~(uint16_t)((ESqrR0 * sigA) >> 16);
        uint16_t recipSqrt16 = r0 + (((uint32_t)r0 * sigma0) >> 25);
        if (!(recipSqrt16 & 0x8000)) {
            recipSqrt16 = 0x8000;
        }
        uint16_t sigZ = ((uint32_t)(sigA << 5) * recipSqrt16) >> 16;
        if (expA) {
            sigZ >>= 1;
        }

        ++sigZ;
        if (!(sigZ & 7)) {
            uint16_t const shiftedSigZ = sigZ >> 1;
            uint16_t const negRem = shiftedSigZ * shiftedSigZ;
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
}
