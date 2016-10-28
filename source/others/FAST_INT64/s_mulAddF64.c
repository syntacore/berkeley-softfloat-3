
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

#include "internals.h"

#include "specialize.h"
#include "softfloat/functions.h"

float64_t
softfloat_mulAddF64(uint64_t uiA, uint64_t uiB, uint64_t uiC, uint8_t op)
{
    if (softfloat_isNaNF64UI(uiA) || softfloat_isNaNF64UI(uiB) || softfloat_isNaNF64UI(uiC)) {
        return ui64_as_f64(softfloat_propagateNaNF64UI(softfloat_propagateNaNF64UI(uiA, uiB), uiC));
    } else {
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

        if (expA == 0x7FF || expB == 0x7FF) {
            /* a or b is inf, product is inf or undefined, check other operand for zero */
            bool const is_product_undefined =
                expA == 0x7FF ? 0 == expB || 0 == sigB :
                0 == expA && 0 == sigA;
            if (!is_product_undefined) {
                /* product is inf */
                uint64_t const uiZ = packToF64UI(signZ, 0x7FF, 0);
                if (expC != 0x7FF) {
                    /* summand c is finite, return product as result */
                    return ui64_as_f64(uiZ);
                }
                /* summand is inf, check for same sign */
                if (signZ == signC) {
                    /* summands are same sign inf */
                    return ui64_as_f64(uiZ);
                }
                /* summands are different sign inf, undefined sum */
            }
            softfloat_raiseFlags(softfloat_flag_invalid);
            return ui64_as_f64(softfloat_propagateNaNF64UI(defaultNaNF64UI, uiC));
        } else if (expC == 0x7FF) {
            return ui64_as_f64(uiC);
        } else {
            struct exp16_sig64 normExpSig;
            int16_t expZ;
            struct uint128 sig128Z;
            uint64_t sigZ;
            int16_t expDiff;
            struct uint128 sig128C;
            int8_t shiftDist;
            if (!expA) {
                if (!sigA) {
                    if (!(expC | sigC) && (signZ != signC)) {
                        return ui64_as_f64(packToF64UI((softfloat_roundingMode == softfloat_round_min), 0, 0));
                    }
                    return ui64_as_f64(uiC);
                }
                normExpSig = softfloat_normSubnormalF64Sig(sigA);
                expA = normExpSig.exp;
                sigA = normExpSig.sig;
            }
            if (!expB) {
                if (!sigB) {
                    if (!(expC | sigC) && (signZ != signC)) {
                        return ui64_as_f64(packToF64UI((softfloat_roundingMode == softfloat_round_min), 0, 0));
                    }
                    return ui64_as_f64(uiC);
                }
                normExpSig = softfloat_normSubnormalF64Sig(sigB);
                expB = normExpSig.exp;
                sigB = normExpSig.sig;
            }

            expZ = expA + expB - 0x3FE;
            sigA = (sigA | UINT64_C(0x0010000000000000)) << 10;
            sigB = (sigB | UINT64_C(0x0010000000000000)) << 10;
            sig128Z = softfloat_mul64To128(sigA, sigB);
            if (sig128Z.v64 < UINT64_C(0x2000000000000000)) {
                --expZ;
                sig128Z =
                    softfloat_add128(
                        sig128Z.v64, sig128Z.v0, sig128Z.v64, sig128Z.v0);
            }
            if (!expC) {
                if (!sigC) {
                    --expZ;
                    return softfloat_roundPackToF64(signZ, expZ, sig128Z.v64 << 1 | (sig128Z.v0 != 0));
                }
                normExpSig = softfloat_normSubnormalF64Sig(sigC);
                expC = normExpSig.exp;
                sigC = normExpSig.sig;
            }
            sigC = (sigC | UINT64_C(0x0010000000000000)) << 9;

            expDiff = expZ - expC;
            if (expDiff < 0) {
                expZ = expC;
                if ((signZ == signC) || (expDiff < -1)) {
                    sig128Z.v64 = softfloat_shiftRightJam64(sig128Z.v64, -expDiff);
                } else {
                    sig128Z =
                        softfloat_shortShiftRightJam128(sig128Z.v64, sig128Z.v0, 1);
                }
            } else if (expDiff) {
                sig128C = softfloat_shiftRightJam128(sigC, 0, expDiff);
            }

            if (signZ == signC) {

                if (expDiff <= 0) {
                    sigZ = (sigC + sig128Z.v64) | (sig128Z.v0 != 0);
                } else {
                    sig128Z = softfloat_add128(sig128Z.v64, sig128Z.v0, sig128C.v64, sig128C.v0);
                    sigZ = sig128Z.v64 | (sig128Z.v0 != 0);
                }
                if (sigZ < UINT64_C(0x4000000000000000)) {
                    --expZ;
                    sigZ <<= 1;
                }
            } else {

                if (expDiff < 0) {
                    signZ = signC;
                    sig128Z = softfloat_sub128(sigC, 0, sig128Z.v64, sig128Z.v0);
                } else if (!expDiff) {
                    sig128Z.v64 = sig128Z.v64 - sigC;
                    if (!(sig128Z.v64 | sig128Z.v0)) {
                        return ui64_as_f64(packToF64UI((softfloat_roundingMode == softfloat_round_min), 0, 0));
                    }
                    if (sig128Z.v64 & UINT64_C(0x8000000000000000)) {
                        signZ = !signZ;
                        sig128Z = softfloat_sub128(0, 0, sig128Z.v64, sig128Z.v0);
                    }
                } else {
                    sig128Z =
                        softfloat_sub128(sig128Z.v64, sig128Z.v0, sig128C.v64, sig128C.v0);
                }

                if (!sig128Z.v64) {
                    expZ -= 64;
                    sig128Z.v64 = sig128Z.v0;
                    sig128Z.v0 = 0;
                }
                shiftDist = softfloat_countLeadingZeros64(sig128Z.v64) - 1;
                expZ -= shiftDist;
                if (shiftDist < 0) {
                    sigZ = softfloat_shortShiftRightJam64(sig128Z.v64, -shiftDist);
                } else {
                    sig128Z =
                        softfloat_shortShiftLeft128(sig128Z.v64, sig128Z.v0, shiftDist);
                    sigZ = sig128Z.v64;
                }
                sigZ |= (sig128Z.v0 != 0);
            }
            return softfloat_roundPackToF64(signZ, expZ, sigZ);
        }
    }
}
