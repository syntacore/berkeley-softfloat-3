
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

#include "target.hpp"

float64_t
f64_sqrt(float64_t const a)
{
    using namespace softfloat::internals;
    uint64_t const uiA = f_as_u_64(a);
    bool const signA = is_sign(uiA);
    int16_t expA = expF64UI(uiA);
    uint64_t sigA = fracF64UI(uiA);

    if (expA == 0x7FF) {
        if (sigA) {
            return u_as_f_64(softfloat_propagateNaNF64UI(uiA, 0));
        }

        if (0 == signA) {
            return a;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_64(defaultNaNF64UI);
    }

    if (signA) {
        if (0 == (expA | sigA)) {
            return a;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_64(defaultNaNF64UI);
    }

    if (0 == expA) {
        if (0 == sigA) {
            return a;
        }

        exp16_sig64 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    /*
    (`sig32Z' is guaranteed to be a lower bound on the square root of
    `sig32A', which makes `sig32Z' also a lower bound on the square root of
    `sigA'.)
    */
    int16_t expZ = ((expA - 0x3FF) >> 1) + 0x3FE;
    expA &= 1;
    sigA |= UINT64_C(0x0010000000000000);
    uint32_t const sig32A = static_cast<uint32_t>(sigA >> 21);
    uint32_t const recipSqrt32 = softfloat_approxRecipSqrt32_1(static_cast<uint32_t>(expA), sig32A);
    uint32_t sig32Z = (static_cast<uint64_t>(sig32A) * recipSqrt32) >> 32;

    if (expA) {
        sigA <<= 8;
        sig32Z >>= 1;
    } else {
        sigA <<= 9;
    }

    uint64_t rem = sigA - static_cast<uint64_t>(sig32Z) * sig32Z;
    uint32_t const q = (static_cast<uint32_t>(rem >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    uint64_t sigZ = (static_cast<uint64_t>(sig32Z) << 32 | 1 << 5) + (static_cast<uint64_t>(q) << 3);

    if ((sigZ & 0x1FF) < 1 << 5) {
        sigZ &= ~static_cast<uint64_t>(0x3F);
        uint64_t const shiftedSigZ = sigZ >> 6;
        rem = (sigA << 52) - shiftedSigZ * shiftedSigZ;

        if (rem & UINT64_C(0x8000000000000000)) {
            --sigZ;
        } else {
            if (rem) {
                sigZ |= 1;
            }
        }
    }

    return softfloat_roundPackToF64(0, expZ, sigZ);
}
