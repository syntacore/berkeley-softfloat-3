
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

float32_t
softfloat_mulAddF32(uint32_t uiA, uint32_t uiB, uint32_t uiC, uint8_t op)
{
    if (softfloat_isNaNF32UI(uiA) || softfloat_isNaNF32UI(uiB) || softfloat_isNaNF32UI(uiC)) {
        return u_as_f_32(softfloat_propagateNaNF32UI(softfloat_propagateNaNF32UI(uiA, uiB), uiC));
    } else {
        bool const signA = signF32UI(uiA);
        int16_t expA = expF32UI(uiA);
        uint32_t sigA = fracF32UI(uiA);
        bool const signB = signF32UI(uiB);
        int16_t expB = expF32UI(uiB);
        uint32_t sigB = fracF32UI(uiB);
        bool const signC = signF32UI(uiC) ^ (op == softfloat_mulAdd_subC);
        int16_t expC = expF32UI(uiC);
        uint32_t sigC = fracF32UI(uiC);
        bool const signProd = signA ^ signB ^ (op == softfloat_mulAdd_subProd);

        if (expA == 0xFF || expB == 0xFF) {
            bool const is_product_undefined = 
                expA == 0xFF ? /* a is inf, check b for zero  */ 0 == expB && 0 == sigB :
                /* expB == 0xFF b is inf, check a for zero*/ 0 == expA && 0 == sigA;
            /* product is inf or undefined */
            if (!is_product_undefined) {
                /* product is inf */
                uint32_t const uiZ = signed_inf_F32UI(signProd);
                if (expC != 0xFF) {
                    /* summand c is finite, return product as result */
                    return u_as_f_32(uiZ);
                }
                /* summand c is inf, check for same sign */
                if (signProd == signC) {
                    /* summands are same sign inf */
                    return u_as_f_32(uiZ);
                }
                /* summands are different sign inf, undefined sum */
            }
            /* product is undefined or undefined sum */
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_32(softfloat_propagateNaNF32UI(defaultNaNF32UI, uiC));
        } else if (expC == 0xFF) {
            return u_as_f_32(uiC);
        } else {
            if (0 == expA) {
                /* a is zero or subnormal */
                if (0 == sigA) {
                    /* a is zero */
                    return 0 == (expC | sigC) && signProd != signC ? signed_zero_F32(softfloat_round_min == softfloat_roundingMode) : u_as_f_32(uiC);
                } else {
                    struct exp16_sig32 const normExpSig = softfloat_normSubnormalF32Sig(sigA);
                    expA = normExpSig.exp;
                    sigA = normExpSig.sig;
                }
            }
            if (0 == expB) {
                /* b is zero or subnormal */
                if (!sigB) {
                    /* b is zero */
                    return 0 == (expC | sigC) && signProd != signC ? signed_zero_F32(softfloat_round_min == softfloat_roundingMode) : u_as_f_32(uiC);
                } else {
                    struct exp16_sig32 const normExpSig = softfloat_normSubnormalF32Sig(sigB);
                    expB = normExpSig.exp;
                    sigB = normExpSig.sig;
                }
            }

            {
                int16_t expProd = expA + expB - 0x7E;
                sigA = (sigA | 0x00800000) << 7;
                sigB = (sigB | 0x00800000) << 7;
                uint64_t sigProd = (uint64_t)sigA * sigB;
                if (sigProd < UINT64_C(0x2000000000000000)) {
                    --expProd;
                    sigProd <<= 1;
                }
                {
                    bool signZ = signProd;
                    int16_t expZ;
                    uint32_t sigZ;
                    if (!expC) {
                        if (!sigC) {
                            sigZ = softfloat_shortShiftRightJam64(sigProd, 31);
                            return softfloat_roundPackToF32(signZ, expProd - 1, sigZ);
                        } else {
                            struct exp16_sig32 const normExpSig = softfloat_normSubnormalF32Sig(sigC);
                            expC = normExpSig.exp;
                            sigC = normExpSig.sig;
                        }
                    }
                    sigC = (sigC | 0x00800000) << 6;

                    {
                        int16_t expDiff = expProd - expC;
                        uint64_t sig64Z;
                        uint64_t sig64C;
                        if (signProd == signC) {
                            if (expDiff <= 0) {
                                expZ = expC;
                                sigZ = sigC + softfloat_shiftRightJam64(sigProd, 32 - expDiff);
                            } else {
                                expZ = expProd;
                                sig64Z = sigProd + softfloat_shiftRightJam64((uint64_t)sigC << 32, expDiff);
                                sigZ = (uint32_t)softfloat_shortShiftRightJam64(sig64Z, 32);
                            }
                            if (sigZ < 0x40000000) {
                                --expZ;
                                sigZ <<= 1;
                            }
                        } else {
                            sig64C = (uint64_t)sigC << 32;
                            if (expDiff < 0) {
                                signZ = signC;
                                expZ = expC;
                                sig64Z = sig64C - softfloat_shiftRightJam64(sigProd, -expDiff);
                            } else if (!expDiff) {
                                expZ = expProd;
                                sig64Z = sigProd - sig64C;
                                if (!sig64Z) {
                                    return signed_zero_F32(softfloat_round_min == softfloat_roundingMode);
                                } else if (sig64Z & INT64_MIN) {
                                    signZ = !signZ;
                                    sig64Z = -(int64_t)sig64Z;
                                }
                            } else {
                                expZ = expProd;
                                sig64Z = sigProd - softfloat_shiftRightJam64(sig64C, expDiff);
                            }
                            {
                                int8_t shiftDist = softfloat_countLeadingZeros64(sig64Z) - 1;
                                expZ -= shiftDist;
                                shiftDist -= 32;
                                sigZ =
                                    shiftDist < 0 ? softfloat_shortShiftRightJam64(sig64Z, -shiftDist) :
                                    (uint32_t)sig64Z << shiftDist;
                            }
                        }
                        return softfloat_roundPackToF32(signZ, expZ, sigZ);
                    }
                }
            }
        }
    }
}

