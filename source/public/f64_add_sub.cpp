
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

using namespace softfloat::internals;

namespace {

static float64_t
addMags(uint64_t const uiA,
        uint64_t const uiB,
        bool const signZ)
{
    using namespace softfloat::internals;
    int16_t const expA = get_exp(uiA);
    uint64_t sigA = get_frac(uiA);
    int16_t const expB = get_exp(uiB);
    uint64_t sigB = get_frac(uiB);
    int16_t const expDiff = expA - expB;

    if (0 == expDiff) {
        if (0 == expA) {
            return u_as_f(uiA + sigB);
        }

        if (0x7FF == expA) {
            return u_as_f(0 != (sigA | sigB) ? propagate_NaN(uiA, uiB) : uiA);
        }

        return round_pack_to_F64(signZ, expA, (UINT64_C(0x0020000000000000) + sigA + sigB) << 9);
    }

    int16_t expZ;
    sigA <<= 9;
    sigB <<= 9;

    if (expDiff < 0) {
        if (expB == 0x7FF) {
            return u_as_f(sigB ? propagate_NaN(uiA, uiB) : pack_to_F64_UI(signZ, 0x7FF, 0));
        }

        expZ = expB;

        if (expA) {
            sigA += UINT64_C(0x2000000000000000);
        } else {
            sigA <<= 1;
        }

        sigA = shift_right_jam_64(sigA, static_cast<uint32_t>(-expDiff));
    } else {
        if (0x7FF == expA) {
            return u_as_f(sigA ? propagate_NaN(uiA, uiB) : uiA);
        }

        expZ = expA;

        if (0 != expB) {
            sigB += UINT64_C(0x2000000000000000);
        } else {
            sigB <<= 1;
        }

        sigB = shift_right_jam_64(sigB, static_cast<uint32_t>(expDiff));
    }

    uint64_t sigZ = UINT64_C(0x2000000000000000) + sigA + sigB;

    if (sigZ < UINT64_C(0x4000000000000000)) {
        --expZ;
        sigZ <<= 1;
    }

    return round_pack_to_F64(signZ, expZ, sigZ);
}

static float64_t
subMags(uint64_t const uiA,
        uint64_t const uiB,
        bool signZ)
{
    int16_t expA = get_exp(uiA);
    uint64_t const sigA = get_frac(uiA);
    int16_t expB = get_exp(uiB);
    uint64_t const sigB = get_frac(uiB);

    int16_t const expDiff = expA - expB;

    if (0 == expDiff) {
        if (0x7FF == expA) {
            if (0 != sigA || 0 != sigB) {
                return u_as_f(propagate_NaN(uiA, uiB));
            }

            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f(defaultNaNF64UI);
        }

        int64_t sigDiff = static_cast<int64_t>(sigA - sigB);

        if (0 == sigDiff) {
            return u_as_f(pack_to_F64_UI(softfloat_round_min == softfloat_get_roundingMode(), 0, 0));
        }

        if (0 != expA) {
            --expA;
        }

        if (sigDiff < 0) {
            signZ = !signZ;
            sigDiff = -sigDiff;
        }

        int8_t const shiftDist = count_leading_zeros(static_cast<uint64_t>(sigDiff)) - 11;
        int16_t const expZ = expA - shiftDist;
        return
            u_as_f(expZ < 0 ?
                      pack_to_F64_UI(signZ, 0, static_cast<uint64_t>(sigDiff << expA)) :
                      pack_to_F64_UI(signZ, expZ, static_cast<uint64_t>(sigDiff << shiftDist)));
    }

    uint64_t const sigA_shifted = sigA << 10;
    uint64_t const sigB_shifted = sigB << 10;

    if (expDiff < 0) {
        signZ = !signZ;

        if (0x7FF == expB) {
            return
                u_as_f(0 != sigB_shifted ?
                          propagate_NaN(uiA, uiB) :
                          pack_to_F64_UI(signZ, 0x7FF, 0));
        }

        return
            norm_round_pack_to_F64(signZ,
                                         expB - 1,
                                         (sigB_shifted | UINT64_C(0x4000000000000000)) -
                                         shift_right_jam_64(sigA_shifted + (expA ? UINT64_C(0x4000000000000000) : sigA_shifted),
                                                 static_cast<uint32_t>(-expDiff)));
    }

    if (0x7FF == expA) {
        return u_as_f(0 == sigA_shifted ? uiA : propagate_NaN(uiA, uiB));
    }

    return
        norm_round_pack_to_F64(signZ,
                                     expA - 1,
                                     (sigA_shifted | UINT64_C(0x4000000000000000)) -
                                     shift_right_jam_64(sigB_shifted + (expB ? UINT64_C(0x4000000000000000) : sigB_shifted), static_cast<uint32_t>(expDiff)));
}

}  // namespace

float64_t
f64_add(float64_t a,
        float64_t b)
{
    using namespace softfloat::internals;
    bool const signA = is_sign(a);
    bool signB = is_sign(b);
    return
        (signA == signB ? addMags : subMags)(f_as_u(a), f_as_u(b), signA);
}

float64_t
f64_sub(float64_t a, float64_t b)
{
    using namespace softfloat::internals;
    bool const signA = is_sign(a);
    bool const signB = is_sign(b);
    return
        (signA == signB ? subMags : addMags)(f_as_u(a), f_as_u(b), signA);
}
