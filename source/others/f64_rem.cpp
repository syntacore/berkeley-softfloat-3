
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

float64_t f64_rem(float64_t a, float64_t b)
{
    uint64_t const uiA = f_as_u_64(a);
    bool const signA = signF64UI(uiA);
    int16_t expA = expF64UI(uiA);
    uint64_t sigA = fracF64UI(uiA);
    uint64_t const uiB = f_as_u_64(b);
    int16_t expB = expF64UI(uiB);
    uint64_t sigB = fracF64UI(uiB);

    if (expA == 0x7FF) {
        if (sigA || ((expB == 0x7FF) && sigB)) {
            return u_as_f_64(softfloat_propagateNaNF64UI(uiA, uiB));
        } else {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_64(defaultNaNF64UI);
        }
    } else if (expB == 0x7FF) {
        return sigB ? u_as_f_64(softfloat_propagateNaNF64UI(uiA, uiB)) : a;
    } else if (expA < expB - 1) {
        return a;
    } else {
        if (!expB) {
            if (!sigB) {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return u_as_f_64(defaultNaNF64UI);
            } else {
                exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(sigB);
                expB = normExpSig.exp;
                sigB = normExpSig.sig;
            }
        }
        if (!expA) {
            if (!sigA) {
                return a;
            } else {
                exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(sigA);
                expA = normExpSig.exp;
                sigA = normExpSig.sig;
            }
        }
        uint64_t rem = sigA | UINT64_C(0x0010000000000000);
        sigB |= UINT64_C(0x0010000000000000);
        int16_t expDiff = expA - expB;
        uint32_t q;
        uint64_t altRem;
        if (expDiff < 1) {
            if (expDiff < -1) {
                return a;
            } else {
                sigB <<= 9;
                if (expDiff) {
                    rem <<= 8;
                    q = 0;
                } else {
                    rem <<= 9;
                    q = (sigB <= rem);
                    if (q) {
                        rem -= sigB;
                    }
                }
            }
        } else {
            uint32_t const recip32 = softfloat_approxRecip32_1(sigB >> 21);
            /*
            Changing the shift of `rem' here requires also changing the initial
            subtraction from `expDiff'.
            */
            rem <<= 9;
            expDiff -= 30;
            /*
            The scale of `sigB' affects how many bits are obtained during each
            cycle of the loop.  Currently this is 29 bits per loop iteration,
            the maximum possible.
            */
            sigB <<= 9;
            uint64_t q64;
            for (;;) {
                q64 = (uint32_t)(rem >> 32) * (uint64_t)recip32;
                if (expDiff < 0) {
                    break;
                }
                q = (q64 + 0x80000000) >> 32;
#ifdef SOFTFLOAT_FAST_INT64
                rem <<= 29;
#else
                rem = (uint64_t)(uint32_t)(rem >> 3) << 32;
#endif
                rem -= q * (uint64_t)sigB;
                if (rem & UINT64_C(0x8000000000000000)) {
                    rem += sigB;
                }
                expDiff -= 29;
            }
            /* `expDiff' cannot be less than -29 here. */
            q = (uint32_t)(q64 >> 32) >> (~expDiff & 31);
            rem = (rem << (expDiff + 30)) - q * (uint64_t)sigB;
            if (rem & UINT64_C(0x8000000000000000)) {
                altRem = rem + sigB;
                goto selectRem;
            }
        }

        do {
            altRem = rem;
            ++q;
            rem -= sigB;
        } while (!(rem & UINT64_C(0x8000000000000000)));
selectRem:
        {
            uint64_t const meanRem = rem + altRem;
            if ((meanRem & UINT64_C(0x8000000000000000)) || (!meanRem && (q & 1))) {
                rem = altRem;
            }
            bool signRem = signA;
            if (rem & UINT64_C(0x8000000000000000)) {
                signRem = !signRem;
                rem = -(int64_t)rem;
            }
            return softfloat_normRoundPackToF64(signRem, expB, rem);
        }
    }
}
