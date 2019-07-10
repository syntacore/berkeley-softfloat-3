
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

namespace {

#ifdef SOFTFLOAT_FAST_DIV64TO32

static inline float32_t
make_result(uint32_t sigA, uint32_t sigB, int16_t expZ, bool const signZ)
{
    using namespace softfloat::internals;
    uint64_t sig64A;
    uint32_t sigZ;

    if (sigA < sigB) {
        --expZ;
        sig64A = static_cast<uint64_t>(sigA) << 31;
    } else {
        sig64A = static_cast<uint64_t>(sigA) << 30;
    }

    sigZ = sig64A / sigB;

    if (!(sigZ & 0x3F)) {
        sigZ |= (static_cast<uint64_t>(sigB) * sigZ != sig64A);
    }

    return softfloat_roundPackToF32(signZ, expZ, sigZ);
}

#else

static inline float32_t
make_result(uint32_t sigA, uint32_t sigB, int16_t expZ, bool const signZ)
{
    using namespace softfloat::internals;
    uint32_t sigZ;
    uint64_t rem;
    if (sigA < sigB) {
        --expZ;
        sigA <<= 8;
    } else {
        sigA <<= 7;
    }

    sigB <<= 8;
    sigZ = (static_cast<uint64_t>(sigA) * softfloat_approxRecip32_1(sigB)) >> 32;

    sigZ += 2;

    if ((sigZ & 0x3F) < 2) {
        sigZ &= ~3;
        rem = (static_cast<uint64_t>(sigA) << 31) - static_cast<uint64_t>(sigZ) * sigB;

        if (rem & UINT64_C(0x8000000000000000)) {
            sigZ -= 4;
        } else {
            if (rem) {
                sigZ |= 1;
            }
        }
    }

    return softfloat_roundPackToF32(signZ, expZ, sigZ);
}

#endif

}  // namespace

float32_t
f32_div(float32_t a,
        float32_t b)
{
    using namespace softfloat::internals;
    uint32_t const uiA = f_as_u_32(a);
    bool const signA = is_sign(uiA);
    int16_t expA = expF32UI(uiA);
    uint32_t sigA = fracF32UI(uiA);
    uint32_t const uiB = f_as_u_32(b);
    bool const signB = is_sign(uiB);
    int16_t expB = expF32UI(uiB);
    uint32_t sigB = fracF32UI(uiB);
    bool const signZ = signA != signB;

    if (0xFF == expA) {
        if (0 != sigA) {
            return u_as_f_32(propagate_NaN(uiA, uiB));
        } else if (0xFF != expB) {
            return signed_inf_F32(signZ);
        } else if (0 != sigB) {
            return u_as_f_32(propagate_NaN(uiA, uiB));
        } else {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_32(defaultNaNF32UI);
        }
    } else if (0xFF == expB) {
        return sigB ? u_as_f_32(propagate_NaN(uiA, uiB)) : signed_zero_F32(signZ);
    } else {
        if (0 == expB) {
            if (0 == sigB) {
                if (0 == (expA | sigA)) {
                    softfloat_raiseFlags(softfloat_flag_invalid);
                    return u_as_f_32(defaultNaNF32UI);
                } else {
                    softfloat_raiseFlags(softfloat_flag_infinite);
                    return signed_inf_F32(signZ);
                }
            } else {
                exp16_sig32 const normExpSig(sigB);
                expB = normExpSig.exp;
                sigB = normExpSig.sig;
            }
        }

        if (0 == expA) {
            if (!sigA) {
                return signed_zero_F32(signZ);
            } else {
                exp16_sig32 const normExpSig(sigA);
                expA = normExpSig.exp;
                sigA = normExpSig.sig;
            }
        }

        return make_result(sigA | 0x00800000, sigB | 0x00800000, expA - expB + 0x7E, signZ);
    }
}
