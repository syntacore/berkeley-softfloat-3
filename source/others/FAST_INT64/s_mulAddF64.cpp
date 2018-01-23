
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

#include "internals.hpp"

#include "target.hpp"
#include "softfloat/functions.h"

namespace softfloat {
namespace internals {

namespace {
static inline float64_t
softfloat_mulAddF64(uint64_t uiA,
                    uint64_t uiB,
                    uint64_t uiC,
                    Mul_add_operations op)
{
    if (softfloat_isNaNF64UI(uiC) || softfloat_isNaNF64UI(uiA) || softfloat_isNaNF64UI(uiB)) {
        return u_as_f_64(softfloat_propagateNaNF64UI(softfloat_propagateNaNF64UI(uiA, uiB), uiC));
    }

    bool const signA = signF64UI(uiA);
    int16_t expA = expF64UI(uiA);
    uint64_t sigA = fracF64UI(uiA);
    bool const signB = signF64UI(uiB);
    int16_t expB = expF64UI(uiB);
    uint64_t sigB = fracF64UI(uiB);
    bool const signC = signF64UI(uiC) ^ (op == softfloat_mulAdd_subC);
    int16_t expC = expF64UI(uiC);
    uint64_t sigC = fracF64UI(uiC);
    bool signZ = signA ^ signB ^ (op == softfloat_mulAdd_subProd);

    if (isInf64UI(uiA) || isInf64UI(uiB)) {
        /* a or b is inf, product is inf or undefined, check other operand for zero */
        bool const is_product_undefined = isInf64UI(uiA) ? isZero64UI(uiB) : isZero64UI(uiA);

        if (!is_product_undefined) {
            /* product is inf */
            uint64_t const uiZ = packToF64UI(signZ, 0x7FF, 0);

            if (expC != 0x7FF) {
                /* summand c is finite, return product as result */
                return u_as_f_64(uiZ);
            }

            /* summand is inf, check for same sign */
            if (signZ == signC) {
                /* summands are same sign inf */
                return u_as_f_64(uiZ);
            }

            /* summands are different sign inf, undefined sum */
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_64(softfloat_propagateNaNF64UI(defaultNaNF64UI, uiC));
    }

    if (0x7FF == expC) {
        // infinity C
        return u_as_f_64(uiC);
    }

    if (0 == expA) {
        // subnormal or zero A
        if (0 == sigA) {
            // zero A
            if (isZero64UI(uiC) && signZ != signC) {
                return u_as_f_64(packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
            }

            return u_as_f_64(uiC);
        }

        exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (0 == expB) {
        if (0 == sigB) {
            if (0 == (expC | sigC) && signZ != signC) {
                return u_as_f_64(packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
            }

            return u_as_f_64(uiC);
        }

        exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int16_t expZ = expA + expB - 0x3FE;
    sigA = (sigA | UINT64_C(0x0010000000000000)) << 10;
    sigB = (sigB | UINT64_C(0x0010000000000000)) << 10;
    uint128 sig128Z = softfloat_mul64To128(sigA, sigB);

    if (sig128Z.v64 < UINT64_C(0x2000000000000000)) {
        --expZ;
        sig128Z = softfloat_add128(sig128Z.v64, sig128Z.v0, sig128Z.v64, sig128Z.v0);
    }

    if (0 == expC) {
        if (0 == sigC) {
            --expZ;
            return softfloat_roundPackToF64(signZ, expZ, sig128Z.v64 << 1 | (sig128Z.v0 != 0));
        }

        exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(sigC);
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
            sig128Z = softfloat_shortShiftRightJam128(sig128Z.v64, sig128Z.v0, 1);
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
            sig128Z = softfloat_add128(sig128Z.v64, sig128Z.v0, sig128C.v64, sig128C.v0);
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
            return u_as_f_64(packToF64UI((softfloat_roundingMode == softfloat_round_min), 0, 0));
        }

        if (sig128Z.v64 & UINT64_C(0x8000000000000000)) {
            signZ = !signZ;
            sig128Z = softfloat_sub128(0, 0, sig128Z.v64, sig128Z.v0);
        }
    } else {
        /* 0 < expDiff */
        assert(0 < expDiff);
        sig128Z = softfloat_sub128(sig128Z.v64, sig128Z.v0, sig128C.v64, sig128C.v0);
    }

    if (0 == sig128Z.v64) {
        expZ -= 64;
        sig128Z.v64 = sig128Z.v0;
        sig128Z.v0 = 0;
    }

    int16_t const shiftDist = softfloat_countLeadingZeros64(sig128Z.v64) - 1;
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

}  // namespace internals
}  // namespace softfloat

float64_t
f64_mulAdd(float64_t a,
           float64_t b,
           float64_t c)
{
    using namespace softfloat::internals;
    return softfloat_mulAddF64(f_as_u_64(a), f_as_u_64(b), f_as_u_64(c), softfloat_mulAdd_madd);
}

