
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
All rights reserved.

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

float32_t
f32_rem(float32_t a,
        float32_t b)
{
    using namespace softfloat;
    uint32_t const uiA = f_as_u_32(a);
    uint32_t const uiB = f_as_u_32(b);

    if (softfloat_isNaNF32UI(uiA) || softfloat_isNaNF32UI(uiB)) {
        return u_as_f_32(softfloat_propagateNaNF32UI(uiA, uiB));
    }

    bool const signA = signF32UI(uiA);
    int16_t expA = expF32UI(uiA);
    uint32_t sigA = fracF32UI(uiA);
    int16_t expB = expF32UI(uiB);
    uint32_t sigB = fracF32UI(uiB);

    if (expA == 0xFF) {
        if (sigA || ((expB == 0xFF) && sigB)) {
            return u_as_f_32(softfloat_propagateNaNF32UI(uiA, uiB));
        }

        /** @todo check */
        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_32(defaultNaNF32UI);
    }

    if (expB == 0xFF) {
        if (sigB) {
            return u_as_f_32(softfloat_propagateNaNF32UI(uiA, uiB));
        }

        /** @todo check */
        return a;
    }

    uint32_t rem;
    int16_t expDiff;
    uint32_t q, recip32, altRem, meanRem;

    if (!expB) {
        if (!sigB) {
            /** @todo check */
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_32(defaultNaNF32UI);
        }

        exp16_sig32 const normExpSig = softfloat_normSubnormalF32Sig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    if (!expA) {
        if (!sigA) {
            return a;
        }

        exp16_sig32 const normExpSig = softfloat_normSubnormalF32Sig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    rem = sigA | 0x00800000;
    sigB |= 0x00800000;
    expDiff = expA - expB;

    if (expDiff < 1) {
        if (expDiff < -1) {
            return a;
        }

        sigB <<= 6;

        if (expDiff) {
            rem <<= 5;
            q = 0;
        } else {
            rem <<= 6;
            q = 0u + !!(sigB <= rem);

            if (0 != q) {
                rem -= sigB;
            }
        }
    } else {
        recip32 = softfloat_approxRecip32_1(sigB << 8);
        /*
         Changing the shift of `rem' here requires also changing the initial
         subtraction from `expDiff'.
        */
        rem <<= 7;
        expDiff -= 31;
        /*
        The scale of `sigB' affects how many bits are obtained during each
        cycle of the loop.  Currently this is 29 bits per loop iteration,
        which is believed to be the maximum possible.
        */
        sigB <<= 6;

        for (;;) {
            q = (rem * static_cast<uint64_t>(recip32)) >> 32;

            if (expDiff < 0) {
                break;
            }

            rem = static_cast<uint32_t>(-static_cast<int32_t>(q * sigB));
            expDiff -= 29;
        }

        /* `expDiff' cannot be less than -30 here. */
        q >>= ~expDiff & 31;
        rem = (rem << (expDiff + 30)) - q * (uint32_t)sigB;
    }

    do {
        altRem = rem;
        ++q;
        rem -= sigB;
    } while (!(rem & 0x80000000));

    meanRem = rem + altRem;

    if ((meanRem & 0x80000000) || (!meanRem && (q & 1))) {
        rem = altRem;
    }

    {
        bool signRem = signA;

        if (0x80000000 <= rem) {
            signRem = !signRem;
            rem = static_cast<uint32_t>(-static_cast<int32_t>(rem));
        }

        return softfloat_normRoundPackToF32(signRem, expB, rem);
    }
}
