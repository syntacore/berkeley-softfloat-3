
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

float16_t
f16_mul(float16_t a,
        float16_t b)
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
        if (sigA || ((expB == 0x1F) && sigB)) {
            return u_as_f_16(softfloat_propagateNaNF16UI(uiA, uiB));
        }

        uint16_t const magBits = static_cast<uint16_t>(expB | sigB);

        if (!magBits) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_16(defaultNaNF16UI);
        }

        return u_as_f_16(packToF16UI(signZ, 0x1F, 0));
    }

    if (expB == 0x1F) {
        if (sigB) {
            return u_as_f_16(softfloat_propagateNaNF16UI(uiA, uiB));
        }

        uint16_t const magBits = static_cast<uint16_t>(expA | sigA);

        if (!magBits) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_16(defaultNaNF16UI);
        }

        return u_as_f_16(packToF16UI(signZ, 0x1F, 0));
    }

    if (!expA) {
        if (!sigA) {
            return u_as_f_16(packToF16UI(signZ, 0, 0));
        }

        exp8_sig16 const normExpSig{sigA};
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (!expB) {
        if (!sigB) {
            return u_as_f_16(packToF16UI(signZ, 0, 0));
        }

        exp8_sig16 const normExpSig{sigB};
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int8_t expZ = expA + expB - 0xF;
    sigA = static_cast<uint16_t>((sigA | 0x0400) << 4);
    sigB = static_cast<uint16_t>((sigB | 0x0400) << 5);
    uint32_t const sig32Z = static_cast<uint32_t>(sigA) * sigB;
    uint16_t sigZ = sig32Z >> 16;

    if (sig32Z & 0xFFFF) {
        sigZ |= 1;
    }

    if (sigZ < 0x4000) {
        --expZ;
        sigZ <<= 1;
    }

    return softfloat_roundPackToF16(signZ, expZ, sigZ);
}
