
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

#include "model.hpp"

namespace {

#if (SOFTFLOAT_FAST_DIV32TO16)

static inline float16_t
make_result(uint16_t const sigA,
            uint16_t const sigB,
            int8_t expZ,
            bool const signZ)
{
    using namespace softfloat::internals;
    uint32_t sig32A;

    if (sigA < sigB) {
        --expZ;
        sig32A = static_cast<uint32_t>(sigA) << 15;
    } else {
        sig32A = static_cast<uint32_t>(sigA) << 14;
    }

    uint16_t sigZ = static_cast<uint16_t>(sig32A / sigB);

    if (0 == (sigZ & 7)) {
        sigZ |= (static_cast<uint32_t>(sigB) * sigZ != sig32A);
    }

    return roundPackToF16(signZ, expZ, sigZ);
}

#else

static inline float16_t
make_result(uint16_t sigA,
            uint16_t const sigB,
            int8_t expZ,
            bool const signZ)
{
    using namespace softfloat::internals;

    if (sigA < sigB) {
        --expZ;
        sigA <<= 5;
    } else {
        sigA <<= 4;
    }

    int const index = sigB >> 6 & 0xF;
    uint16_t const r0 =
        approxRecip_1k0s[index] -
        ((static_cast<uint32_t>(approxRecip_1k1s[index]) * (sigB & 0x3F)) >> 10);
    uint16_t sigZ = (static_cast<uint32_t>(sigA) * r0) >> 16;
    uint16_t rem = (sigA << 10) - sigZ * sigB;
    sigZ += (rem * static_cast<uint32_t>(r0)) >> 26;

    ++sigZ;

    if (0 == (sigZ & 7)) {
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

    return roundPackToF16(signZ, expZ, sigZ);
}

#endif

}  // namespace

float16_t
f16_div(float16_t const a,
        float16_t const b)
{
    using namespace softfloat::internals;
    bool const signA = is_sign(a);
    int8_t expA = get_exp(a);
    uint16_t sigA = get_frac(a);
    bool const signB = is_sign(b);
    int8_t expB = get_exp(b);
    uint16_t sigB = get_frac(b);
    bool const signZ = signA != signB;

    if (is_NaN(a) || is_NaN(b)) {
        return propagate_NaN(a, b);
    }

    if (is_inf(a)) {
        if (is_finite(b)) {
            return u_as_f(packToF16UI(signZ, 0x1F, 0));
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f(defaultNaNF16UI);
    }

    if (is_inf(b)) {
        return u_as_f(packToF16UI(signZ, 0, 0));
    }

    if (is_zero(b)) {
        if (is_zero(a)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f(defaultNaNF16UI);
        }

        softfloat_raiseFlags(softfloat_flag_infinite);
        return u_as_f(packToF16UI(signZ, 0x1F, 0));
    }

    if (is_zero(a)) {
        return u_as_f(packToF16UI(signZ, 0, 0));
    }

    if (0 == expB) {
        exp8_sig16 const normExpSig{sigB};
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    if (0 == expA) {
        exp8_sig16 const normExpSig{sigA};
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    return
        make_result(static_cast<uint16_t>(UINT16_C(0x0400) | sigA),
                    static_cast<uint16_t>(UINT16_C(0x0400) | sigB),
                    static_cast<int8_t>(expA - expB + 0xEu),
                    signZ);
}
