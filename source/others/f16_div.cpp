
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
#include "target.hpp"

float16_t
f16_div(float16_t const a,
        float16_t const b)
{
    using namespace softfloat::internals;
    uint16_t const uiA = f_as_u_16(a);
    bool const signA = signF16UI(uiA);
    int8_t expA = expF16UI(uiA);
    uint16_t sigA = fracF16UI(uiA);
    uint16_t const uiB = f_as_u_16(b);
    bool const signB = signF16UI(uiB);
    int8_t expB = expF16UI(uiB);
    uint16_t sigB = fracF16UI(uiB);
    bool const signZ = signA ^ signB;

    if (expA == 0x1F) {
        if (sigA) {
            return u_as_f_16(softfloat_propagateNaNF16UI(uiA, uiB));
        }

        if (expB == 0x1F) {
            if (sigB) {
                return u_as_f_16(softfloat_propagateNaNF16UI(uiA, uiB));
            }

            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_16(defaultNaNF16UI);
        }

        return u_as_f_16(packToF16UI(signZ, 0x1F, 0));
    }

    if (expB == 0x1F) {
        if (sigB) {
            return u_as_f_16(softfloat_propagateNaNF16UI(uiA, uiB));
        }

        return u_as_f_16(packToF16UI(signZ, 0, 0));
    }

    if (!expB) {
        if (!sigB) {
            if (!(expA | sigA)) {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return u_as_f_16(defaultNaNF16UI);
            }

            softfloat_raiseFlags(softfloat_flag_infinite);
            return u_as_f_16(packToF16UI(signZ, 0x1F, 0));
        }

        exp8_sig16 const normExpSig = softfloat_normSubnormalF16Sig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    if (!expA) {
        if (!sigA) {
            return u_as_f_16(packToF16UI(signZ, 0, 0));
        }

        exp8_sig16 const normExpSig = softfloat_normSubnormalF16Sig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int8_t expZ = expA - expB + 0xE;
    sigA |= 0x0400u;
    sigB |= 0x0400u;
#ifdef SOFTFLOAT_FAST_DIV32TO16
    uint32_t sig32A;

    if (sigA < sigB) {
        --expZ;
        sig32A = static_cast<uint32_t>(sigA) << 15;
    } else {
        sig32A = static_cast<uint32_t>(sigA) << 14;
    }

    uint16_t sigZ = static_cast<uint16_t>(sig32A / sigB);

    if (!(sigZ & 7)) {
        sigZ |= (static_cast<uint32_t>(sigB) * sigZ != sig32A);
    }

#else

    if (sigA < sigB) {
        --expZ;
        sigA <<= 5;
    } else {
        sigA <<= 4;
    }

    int const index = sigB >> 6 & 0xF;
    uint16_t const r0 =
        softfloat_approxRecip_1k0s[index] - 
        ((static_cast<uint32_t>(softfloat_approxRecip_1k1s[index]) * (sigB & 0x3F)) >> 10);
    uint16_t sigZ = (static_cast<uint32_t>(sigA) * r0) >> 16;
    uint16_t rem = (sigA << 10) - sigZ * sigB;
    sigZ += (rem * static_cast<uint32_t>(r0)) >> 26;

    ++sigZ;

    if (!(sigZ & 7)) {
        sigZ &= ~1;
        rem = (sigA << 10) - sigZ * sigB;

        if (rem & 0x8000) {
            sigZ -= 2;
        } else {
            if (rem) {
                sigZ |= 1;
            }
        }
    }

#endif
    return softfloat_roundPackToF16(signZ, expZ, sigZ);
}
