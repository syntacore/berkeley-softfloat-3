
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

#include "model.hpp"

#if (SOFTFLOAT_FAST_INT64)

namespace {

#if (SOFTFLOAT_FAST_DIV64TO32)

static inline float32_t
make_result(uint32_t sigA, uint32_t sigB, int16_t expZ, bool const signZ)
{
    using namespace softfloat::internals;

    uint64_t sig64A;

    if (sigA < sigB) {
        --expZ;
        sig64A = static_cast<uint64_t>(sigA) << 31;
    } else {
        sig64A = static_cast<uint64_t>(sigA) << 30;
    }

    uint32_t sigZ = static_cast<uint32_t>(sig64A / sigB);

    if (0 == (sigZ & 0x3F)) {
        sigZ |= !!(static_cast<uint64_t>(sigB) * sigZ != sig64A);
    }

    return round_pack_to_F32(signZ, expZ, sigZ);
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
    sigZ = (static_cast<uint64_t>(sigA) * approx_recip_32_1(sigB)) >> 32;

    sigZ += 2;

    if ((sigZ & 0x3F) < 2) {
        sigZ &= ~3;
        rem = (static_cast<uint64_t>(sigA) << 31) - static_cast<uint64_t>(sigZ) * sigB;

        if (rem & UINT64_C(0x8000000000000000)) {
            sigZ -= 4;
        } else if (rem) {
            sigZ |= 1;
        }
    }

    return round_pack_to_F32(signZ, expZ, sigZ);
}

#endif

}  // namespace

float32_t
f32_div(float32_t const a,
        float32_t const b)
{
    using namespace softfloat::internals;
    bool const signA = is_sign(a);
    int16_t expA = get_exp(a);
    uint32_t sigA = get_frac(a);
    bool const signB = is_sign(b);
    int16_t expB = get_exp(b);
    uint32_t sigB = get_frac(b);
    bool const signZ = signA != signB;

    if (is_NaN(a) || is_NaN(b)) {
        return propagate_NaN(a, b);
    }

    if (is_inf(a)) {
        if (is_finite(b)) {
            /**
            @todo check inf / 0 case
            */
            return make_signed_inf<float32_t>(signZ);
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f(defaultNaNF32UI);
    }

    if (is_inf(b)) {
        return make_signed_zero<float32_t>(signZ);
    }

    if (is_zero(b)) {
        if (is_zero(a)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f(defaultNaNF32UI);
        }

        softfloat_raiseFlags(softfloat_flag_infinite);
        return make_signed_inf<float32_t>(signZ);
    }

    if (is_zero(a)) {
        return make_signed_zero<float32_t>(signZ);
    }

    if (0 == expB) {
        exp16_sig32 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    if (0 == expA) {
        exp16_sig32 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    return make_result(sigA | 0x00800000, sigB | 0x00800000, expA - expB + 0x7E, signZ);
}

#else

float32_t
f32_div(float32_t const a,
        float32_t const b)
{
    using namespace softfloat::internals::slow_int64;

    if (is_NaN(a) || is_NaN(b)) {
        return propagate_NaN(a, b);
    }

    bool const signA = is_sign(a);
    bool const signB = is_sign(b);
    bool const signZ = signA != signB;

    if (is_inf(a)) {
        if (is_inf(b)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f(defaultNaNF32UI);
        }

        return make_signed_inf<float32_t>(signZ);
    }

    if (is_inf(b)) {
        return  make_signed_zero<float32_t>(signZ);
    }

    if (is_zero(b)) {
        if (is_zero(a)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f(defaultNaNF32UI);
        }

        softfloat_raiseFlags(softfloat_flag_infinite);
        return make_signed_inf<float32_t>(signZ);
    }

    if (is_zero(a)) {
        return make_signed_zero<float32_t>(signZ);
    }

    int16_t expB = get_exp(b);
    uint32_t sigB = get_frac(b);

    if (0 == expB) {
        exp16_sig32 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int16_t expA = get_exp(a);
    uint32_t sigA = get_frac(a);

    if (0 == expA) {
        exp16_sig32 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    {
        int16_t expZ = expA - expB + 0x7E;
        sigA |= 0x00800000;
        sigB |= 0x00800000;

#if (SOFTFLOAT_FAST_DIV64TO32)

        uint64_t sig64A;

        if (sigA < sigB) {
            --expZ;
            sig64A = static_cast<uint64_t>(sigA) << 31;
        } else {
            sig64A = static_cast<uint64_t>(sigA) << 30;
        }

        uint32_t sigZ = sig64A / sigB;

        if (0 == (sigZ & 0x3F)) {
            sigZ |= (static_cast<uint64_t>(sigB) * sigZ != sig64A);
        }

#else

        if (sigA < sigB) {
            --expZ;
            sigA <<= 8;
        } else {
            sigA <<= 7;
        }

        sigB <<= 8;
        uint32_t sigZ = (static_cast<uint64_t>(sigA) * approx_recip_32_1(sigB)) >> 32;

        sigZ += 2;

        if ((sigZ & 0x3F) < 2) {
            sigZ &= ~3;
            uint64_t const rem = (static_cast<uint64_t>(sigA) << 32) - static_cast<uint64_t>(sigZ << 1) * sigB;

            if (0 != (rem & UINT64_C(0x8000000000000000))) {
                sigZ -= 4;
            } else {
                if (0 != rem) {
                    sigZ |= 1;
                }
            }
        }

#endif
        return round_pack_to_F32(signZ, expZ, sigZ);
    }
}

#endif
