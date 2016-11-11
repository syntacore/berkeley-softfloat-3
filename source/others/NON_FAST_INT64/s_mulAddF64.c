
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
    if (softfloat_isNaNF64UI(uiC) || softfloat_isNaNF64UI(uiA) || softfloat_isNaNF64UI(uiB)) {
        return u_as_f_64(softfloat_propagateNaNF64UI(softfloat_propagateNaNF64UI(uiA, uiB), uiC));
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

        if (isInf64UI(uiA) || isInf64UI(uiB)) {
            /* a or b is inf, product is inf or undefined, check other operand for zero */
            bool const is_product_undefined =
                isInf64UI(uiA) ? isZero64UI(uiB) : isZero64UI(uiA);
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
        } else if (expC == 0x7FF) {
            return u_as_f_64(uiC);
        } else {
            uint32_t sig128C[4];
            int16_t expZ;
            uint32_t sig128Z[4];
            uint64_t sigZ;
            int16_t shiftDist;
            int16_t expDiff;
            if (!expA) {
                if (!sigA) {
                    return u_as_f_64(0 == (expC | sigC) && signZ != signC ? packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0) : uiC);
                }
                struct exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(sigA);
                expA = normExpSig.exp;
                sigA = normExpSig.sig;
            }
            if (!expB) {
                if (!sigB) {
                    return u_as_f_64(0 == (expC | sigC) && signZ != signC ? packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0) : uiC);
                }
                struct exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(sigB);
                expB = normExpSig.exp;
                sigB = normExpSig.sig;
            }

            expZ = expA + expB - 0x3FE;
            sigA = (sigA | UINT64_C(0x0010000000000000)) << 10;
            sigB = (sigB | UINT64_C(0x0010000000000000)) << 11;
            softfloat_mul64To128M(sigA, sigB, sig128Z);
            sigZ =
                (uint64_t)sig128Z[indexWord(4, 3)] << 32 | sig128Z[indexWord(4, 2)];
            shiftDist = 0;
            if (!(sigZ & UINT64_C(0x4000000000000000))) {
                --expZ;
                shiftDist = -1;
            }
            if (!expC) {
                if (!sigC) {
                    if (shiftDist) {
                        sigZ <<= 1;
                    }
                    return
                        softfloat_roundPackToF64(signZ, expZ - 1,
                                                 sigZ | (sig128Z[indexWord(4, 1)] || sig128Z[indexWord(4, 0)]));
                }
                struct exp16_sig64 const normExpSig = softfloat_normSubnormalF64Sig(sigC);
                expC = normExpSig.exp;
                sigC = normExpSig.sig;
            }
            sigC = (sigC | UINT64_C(0x0010000000000000)) << 10;

            expDiff = expZ - expC;
            if (expDiff < 0) {
                expZ = expC;
                if ((signZ == signC) || (expDiff < -1)) {
                    shiftDist -= expDiff;
                    if (shiftDist) {
                        sigZ = softfloat_shiftRightJam64(sigZ, shiftDist);
                    }
                } else {
                    if (!shiftDist) {
                        softfloat_shortShiftRight128M(sig128Z, 1, sig128Z);
                    }
                }
            } else {
                if (shiftDist) {
                    softfloat_add128M(sig128Z, sig128Z, sig128Z);
                }
                if (!expDiff) {
                    sigZ = (uint64_t)sig128Z[indexWord(4, 3)] << 32 | sig128Z[indexWord(4, 2)];
                } else {
                    sig128C[indexWord(4, 3)] = sigC >> 32;
                    sig128C[indexWord(4, 2)] = (uint32_t)sigC;
                    sig128C[indexWord(4, 1)] = 0;
                    sig128C[indexWord(4, 0)] = 0;
                    softfloat_shiftRightJam128M(sig128C, expDiff, sig128C);
                }
            }

            if (signZ == signC) {
                if (expDiff <= 0) {
                    sigZ += sigC;
                } else {
                    softfloat_add128M(sig128Z, sig128C, sig128Z);
                    sigZ =
                        (uint64_t)sig128Z[indexWord(4, 3)] << 32 |
                        sig128Z[indexWord(4, 2)];
                }
                if (sigZ & UINT64_C(0x8000000000000000)) {
                    ++expZ;
                    sigZ = softfloat_shortShiftRightJam64(sigZ, 1);
                }
            } else {
                if (expDiff < 0) {
                    signZ = signC;
                    if (expDiff < -1) {
                        sigZ = sigC - sigZ;
                        if (sig128Z[indexWord(4, 1)] || sig128Z[indexWord(4, 0)]) {
                            sigZ = (sigZ - 1) | 1;
                        }
                        if (!(sigZ & UINT64_C(0x4000000000000000))) {
                            --expZ;
                            sigZ <<= 1;
                        }
                        return softfloat_roundPackToF64(signZ, expZ - 1, sigZ);
                    } else {
                        sig128C[indexWord(4, 3)] = sigC >> 32;
                        sig128C[indexWord(4, 2)] = (uint32_t)sigC;
                        sig128C[indexWord(4, 1)] = 0;
                        sig128C[indexWord(4, 0)] = 0;
                        softfloat_sub128M(sig128C, sig128Z, sig128Z);
                    }
                } else if (!expDiff) {
                    sigZ -= sigC;
                    if (0 == sigZ && !sig128Z[indexWord(4, 1)] && !sig128Z[indexWord(4, 0)]) {
                        return u_as_f_64(packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
                    }
                    sig128Z[indexWord(4, 3)] = sigZ >> 32;
                    sig128Z[indexWord(4, 2)] = (uint32_t)sigZ;
                    if (sigZ & INT64_MIN) {
                        signZ = !signZ;
                        softfloat_negX128M(sig128Z);
                    }
                } else {
                    softfloat_sub128M(sig128Z, sig128C, sig128Z);
                    if (1 < expDiff) {
                        sigZ = (uint64_t)sig128Z[indexWord(4, 3)] << 32 | sig128Z[indexWord(4, 2)];
                        if (!(sigZ & UINT64_C(0x4000000000000000))) {
                            --expZ;
                            sigZ <<= 1;
                        }
                        return
                            softfloat_roundPackToF64(signZ, expZ - 1,
                                                     sigZ | (sig128Z[indexWord(4, 1)] || sig128Z[indexWord(4, 0)]));
                    }
                }

                shiftDist = 0;
                sigZ = (uint64_t)sig128Z[indexWord(4, 3)] << 32 | sig128Z[indexWord(4, 2)];
                if (!sigZ) {
                    shiftDist = 64;
                    sigZ = (uint64_t)sig128Z[indexWord(4, 1)] << 32 | sig128Z[indexWord(4, 0)];
                }
                shiftDist += softfloat_countLeadingZeros64(sigZ) - 1;
                if (shiftDist) {
                    expZ -= shiftDist;
                    softfloat_shiftLeft128M(sig128Z, shiftDist, sig128Z);
                    sigZ = (uint64_t)sig128Z[indexWord(4, 3)] << 32 | sig128Z[indexWord(4, 2)];
                }
            }
            return
                softfloat_roundPackToF64(signZ, expZ - 1,
                                         sigZ | (sig128Z[indexWord(4, 1)] || sig128Z[indexWord(4, 0)]));
        }
    }
}
