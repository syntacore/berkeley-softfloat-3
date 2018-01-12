
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

#include <cassert>

float16_t
f16_rem(float16_t a, float16_t b)
{
    uint16_t const uiA = f_as_u_16(a);
    bool const signA = signF16UI(uiA);
    int8_t expA = expF16UI(uiA);
    uint16_t sigA = fracF16UI(uiA);
    uint16_t const uiB = f_as_u_16(b);
    int8_t expB = expF16UI(uiB);
    uint16_t sigB = fracF16UI(uiB);

    if (expA == 0x1F) {
        if (sigA || ((expB == 0x1F) && sigB)) {
            return u_as_f_16(softfloat_propagateNaNF16UI(uiA, uiB));
        } else {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_16(defaultNaNF16UI);
        }
    } else if (expB == 0x1F) {
        return sigB ? u_as_f_16(softfloat_propagateNaNF16UI(uiA, uiB)) : a;
    } else {
        if (!expB) {
            if (!sigB) {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return u_as_f_16(defaultNaNF16UI);
            } else {
                struct exp8_sig16 const normExpSig = softfloat_normSubnormalF16Sig(sigB);
                expB = normExpSig.exp;
                sigB = normExpSig.sig;
            }
        }
        if (!expA) {
            if (!sigA) {
                return a;
            } else {
                struct exp8_sig16 const normExpSig = softfloat_normSubnormalF16Sig(sigA);
                expA = normExpSig.exp;
                sigA = normExpSig.sig;
            }
        }

        uint16_t rem = sigA | 0x0400;
        sigB |= 0x0400;
        int8_t expDiff = expA - expB;
        uint16_t q;
        if (expDiff < 1) {
            if (expDiff < -1) {
                return a;
            } else {
                sigB <<= 3;
                if (expDiff) {
                    rem <<= 2;
                    q = 0;
                } else {
                    rem <<= 3;
                    q = (sigB <= rem);
                    if (q) {
                        rem -= sigB;
                    }
                }
            }
        } else {
            uint32_t const recip32 = softfloat_approxRecip32_1((uint32_t)sigB << 21);
            /*
            Changing the shift of `rem' here requires also changing the initial
            subtraction from `expDiff'.
            */
            rem <<= 4;
            expDiff -= 31;
            /*
            The scale of `sigB' affects how many bits are obtained during each
            cycle of the loop.  Currently this is 29 bits per loop iteration,
            which is believed to be the maximum possible.
            */
            sigB <<= 3;
            uint32_t q32;
            for (;;) {
                q32 = (rem * (uint64_t)recip32) >> 16;
                if (expDiff < 0) {
                    break;
                }
                rem = -(int16_t)((uint16_t)q32 * sigB);
                expDiff -= 29;
            }
            /*`expDiff' cannot be less than -30 here.*/
            assert(-30 <= expDiff);
            q32 >>= ~expDiff & 31;
            /** @todo Warning	C4242	'=': conversion from 'uint32_t' to 'uint16_t', possible loss of data*/
            q = q32;
            rem = (rem << (expDiff + 30)) - q * sigB;
        }

        uint16_t altRem;
        do {
            altRem = rem;
            ++q;
            rem -= sigB;
        } while (!(rem & 0x8000));
        uint16_t const meanRem = rem + altRem;
        if ((meanRem & 0x8000) || (!meanRem && (q & 1))) {
            rem = altRem;
        }
        bool signRem = signA;
        if (0x8000 <= rem) {
            signRem = !signRem;
            rem = -rem;
        }
        return softfloat_normRoundPackToF16(signRem, expB, rem);
    }
}
