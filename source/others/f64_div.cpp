
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

float64_t
f64_div(float64_t a, float64_t b)
{
    uint64_t const uiA = f_as_u_64(a);
    bool const signA = signF64UI(uiA);
    int16_t expA = expF64UI(uiA);
    uint64_t sigA = fracF64UI(uiA);
    uint64_t uiB = f_as_u_64(b);
    bool const signB = signF64UI(uiB);
    int16_t expB = expF64UI(uiB);
    uint64_t sigB = fracF64UI(uiB);
    bool const signZ = signA ^ signB;

    if (expA == 0x7FF) {
        if (sigA) {
            return u_as_f_64(softfloat_propagateNaNF64UI(uiA, uiB));
        } else if (expB == 0x7FF) {
            if (sigB) {
                return u_as_f_64(softfloat_propagateNaNF64UI(uiA, uiB));
            } else {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return u_as_f_64(defaultNaNF64UI);
            }
        } else {
            return u_as_f_64(packToF64UI(signZ, 0x7FF, 0));
        }
    } else if (expB == 0x7FF) {
        return u_as_f_64(sigB ? softfloat_propagateNaNF64UI(uiA, uiB) : packToF64UI(signZ, 0, 0));
    } else {
        if (!expB) {
            if (!sigB) {
                if (!(expA | sigA)) {
                    softfloat_raiseFlags(softfloat_flag_invalid);
                    return u_as_f_64(defaultNaNF64UI);
                } else {
                    softfloat_raiseFlags(softfloat_flag_infinite);
                    return u_as_f_64(packToF64UI(signZ, 0x7FF, 0));
                }
            } else {
                struct exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(sigB);
                expB = normExpSig.exp;
                sigB = normExpSig.sig;
            }
        }
        if (!expA) {
            if (!sigA) {
                return u_as_f_64(packToF64UI(signZ, 0, 0));
            } else {
                struct exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(sigA);
                expA = normExpSig.exp;
                sigA = normExpSig.sig;
            }
        }

        int16_t expZ = expA - expB + 0x3FE;
        sigA |= UINT64_C(0x0010000000000000);
        sigB |= UINT64_C(0x0010000000000000);
        if (sigA < sigB) {
            --expZ;
            sigA <<= 11;
        } else {
            sigA <<= 10;
        }
        sigB <<= 11;
        uint32_t const recip32 = softfloat_approxRecip32_1(sigB >> 32) - 2;
        uint32_t const sig32Z = ((uint32_t)(sigA >> 32) * (uint64_t)recip32) >> 32;
        uint32_t doubleTerm = sig32Z << 1;
        uint64_t rem =
            ((sigA - (uint64_t)doubleTerm * (uint32_t)(sigB >> 32)) << 28) -
            (uint64_t)doubleTerm * ((uint32_t)sigB >> 4);
        uint32_t q = (((uint32_t)(rem >> 32) * (uint64_t)recip32) >> 32) + 4;
        uint64_t sigZ = ((uint64_t)sig32Z << 32) + ((uint64_t)q << 4);

        if ((sigZ & 0x1FF) < 4 << 4) {
            q &= ~7;
            sigZ &= ~(uint64_t)0x7F;
            doubleTerm = q << 1;
            rem =
                ((rem - (uint64_t)doubleTerm * (uint32_t)(sigB >> 32)) << 28)
                - (uint64_t)doubleTerm * ((uint32_t)sigB >> 4);
            if (rem & UINT64_C(0x8000000000000000)) {
                sigZ -= 1 << 7;
            } else {
                if (rem) {
                    sigZ |= 1;
                }
            }
        }
        return softfloat_roundPackToF64(signZ, expZ, sigZ);
    }
}
