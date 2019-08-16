
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

#ifndef SOFTFLOAT_FAST_INT64
#error Fast int64_t operations only
#endif

using namespace softfloat::internals;

namespace {

static inline float64_t
mulAdd(uint64_t const& uiA,
       uint64_t const& uiB,
       uint64_t const& uiC)
{
    if (is_sNaN(uiC) || is_NaN(uiA) || is_NaN(uiB)) {
        return to_float(propagate_NaN(propagate_NaN(uiA, uiB), uiC));
    }

    bool const signA = is_sign(uiA);
    bool const signB = is_sign(uiB);
    bool const signC = is_sign(uiC);
    int16_t expC = get_exp(uiC);
    bool signZ = signA != signB;

    if (is_inf(uiA) || is_inf(uiB)) {
        /* a or b is inf, product is inf or undefined, check other operand for zero */
        bool const is_product_undefined = is_inf(uiA) ? is_zero(uiB) : is_zero(uiA);

        if (is_product_undefined) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return to_float(propagate_NaN(defaultNaNF64UI, uiC));
        }

        /* product is inf */
        uint64_t const uiZ = packToF64UI(signZ, 0x7FF, 0);

        if (expC != 0x7FF) {
            /* summand c is finite, return product as result */
            return to_float(uiZ);
        }

        /* if summand is inf, check for same sign */
        if (!is_NaN(uiC) && signZ == signC) {
            /* summands are same sign inf */
            return to_float(uiZ);
        }

        /* summands are different sign inf or NaN, undefined sum */
        softfloat_raiseFlags(softfloat_flag_invalid);
        return to_float(propagate_NaN(defaultNaNF64UI, uiC));
    }

    if (is_NaN(uiC)) {
        return to_float(propagate_NaN(defaultNaNF64UI, uiC));
    }

    if (0x7FF == expC) {
        // infinity C
        return to_float(uiC);
    }

    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
    int16_t expA = get_exp(uiA);
    uint64_t sigA = get_frac(uiA);

    if (0 == expA) {
        // subnormal or zero A
        if (0 == sigA) {
            // zero A
            if (is_zero(uiC) && signZ != signC) {
                return to_float(packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
            }

            return to_float(uiC);
        }

        exp16_sig64 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int16_t expB = get_exp(uiB);
    uint64_t sigB = get_frac(uiB);
    uint64_t sigC = get_frac(uiC);

    if (0 == expB) {
        if (0 == sigB) {
            if (0 == (expC | sigC) && signZ != signC) {
                return to_float(packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
            }

            return to_float(uiC);
        }

        exp16_sig64 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int16_t expZ = expA + expB - 0x3FE;
    sigA = (sigA | UINT64_C(0x0010000000000000)) << 10;
    sigB = (sigB | UINT64_C(0x0010000000000000)) << 10;
    uint128 sig128Z = softfloat_mul64To128(sigA, sigB);

    if (sig128Z.v64 < UINT64_C(0x2000000000000000)) {
        --expZ;
        sig128Z = softfloat_add128(sig128Z, sig128Z);
    }

    if (0 == expC) {
        if (0 == sigC) {
            --expZ;
            return softfloat_roundPackToF64(signZ, expZ, sig128Z.v64 << 1 | (sig128Z.v0 != 0));
        }

        exp16_sig64 const normExpSig(sigC);
        expC = normExpSig.exp;
        sigC = normExpSig.sig;
    }

    sigC = (sigC | UINT64_C(0x0010000000000000)) << 9;

    int16_t const expDiff = expZ - expC;

    uint128 sig128C;

    if (expDiff < 0) {
        expZ = expC;

        if (signZ == signC || expDiff < -1) {
            sig128Z.v64 = softfloat_shiftRightJam64(sig128Z.v64, static_cast<uint32_t>(-expDiff));
        } else {
            sig128Z = softfloat_shortShiftRightJam128(sig128Z, 1);
        }
    } else if (0 < expDiff) {
        sig128C = softfloat_shiftRightJam128(sigC, 0u, static_cast<uint32_t>(expDiff));
    }

    if (signZ == signC) {
        uint64_t sigZ;

        if (expDiff <= 0) {
            sigZ = (sigC + sig128Z.v64) | (sig128Z.v0 != 0);
        } else {
            /* 0 < expDiff */
            assert(0 < expDiff);
            sig128Z = softfloat_add128(sig128Z, sig128C);
            sigZ = sig128Z.v64 | (sig128Z.v0 != 0);
        }

        if (sigZ < UINT64_C(0x4000000000000000)) {
            --expZ;
            sigZ <<= 1;
        }

        return softfloat_roundPackToF64(signZ, expZ, sigZ);
    }

    if (expDiff < 0) {
        signZ = signC;
        sig128Z = softfloat_sub128(sigC, 0, sig128Z.v64, sig128Z.v0);
    } else if (0 == expDiff) {
        sig128Z.v64 = sig128Z.v64 - sigC;

        if (0 == (sig128Z.v64 | sig128Z.v0)) {
            return to_float(packToF64UI(softfloat_round_min == softfloat_roundingMode, 0, 0));
        }

        if (sig128Z.v64 & UINT64_C(0x8000000000000000)) {
            signZ = !signZ;
            sig128Z = softfloat_sub128(0, 0, sig128Z.v64, sig128Z.v0);
        }
    } else {
        /* 0 < expDiff */
        assert(0 < expDiff);
        sig128Z = softfloat_sub128(sig128Z, sig128C);
    }

    if (0 == sig128Z.v64) {
        expZ -= 64;
        sig128Z.v64 = sig128Z.v0;
        sig128Z.v0 = 0;
    }

    int16_t const shiftDist = count_leading_zeros(sig128Z.v64) - 1;
    expZ -= shiftDist;

    if (shiftDist < 0) {
        uint64_t const sigZ = softfloat_shortShiftRightJam64(sig128Z.v64, static_cast<uint8_t>(-shiftDist));
        return
            softfloat_roundPackToF64(signZ,
                                     expZ,
                                     sigZ | !!(0 != sig128Z.v0));
    }

    uint128 const sig128Z_1 = softfloat_shortShiftLeft128(sig128Z.v64, sig128Z.v0, static_cast<uint8_t>(shiftDist));
    return softfloat_roundPackToF64(signZ,
                                    expZ,
                                    sig128Z_1.v64 | !!(0 != sig128Z_1.v0));
}
}

float64_t
f64_mulAdd(float64_t a,
           float64_t b,
           float64_t c)
{
    using namespace softfloat::internals;
    return mulAdd(f_as_u(a), f_as_u(b), f_as_u(c));
}

