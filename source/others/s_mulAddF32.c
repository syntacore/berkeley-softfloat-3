
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

static __inline float32_t ui32_as_f32(uint32_t v)
{
    return *(float32_t const*)&v;
}

float32_t
softfloat_mulAddF32(uint32_t uiA, uint32_t uiB, uint32_t uiC, uint8_t op)
{
    bool const signA = signF32UI(uiA);
    int16_t expA = expF32UI(uiA);
    uint32_t sigA = fracF32UI(uiA);
    bool const signB = signF32UI(uiB);
    int16_t expB = expF32UI(uiB);
    uint32_t sigB = fracF32UI(uiB);
    bool signC = signF32UI(uiC) ^ (op == softfloat_mulAdd_subC);
    int16_t expC = expF32UI(uiC);
    uint32_t sigC = fracF32UI(uiC);
    bool const signProd = signA ^ signB ^ (op == softfloat_mulAdd_subProd);

    if (expA == 0xFF) {
        if (sigA || (expB == 0xFF && sigB)) {
            return ui32_as_f32(softfloat_propagateNaNF32UI(softfloat_propagateNaNF32UI(uiA, uiB), uiC));
        } else {
            uint32_t const magBits = expB | sigB;
            if (magBits) {
                uint32_t const uiZ = packToF32UI(signProd, 0xFF, 0);
                if (expC != 0xFF) {
                    return ui32_as_f32(uiZ);
                } else if (sigC) {
                    return ui32_as_f32(softfloat_propagateNaNF32UI(uiZ, uiC));
                } else if (signProd == signC) {
                    return ui32_as_f32(uiZ);
                }
            }
            softfloat_raiseFlags(softfloat_flag_invalid);
            return ui32_as_f32(softfloat_propagateNaNF32UI(defaultNaNF32UI, uiC));
        }
    }
    if (expB == 0xFF) {
        if (sigB) {
            return ui32_as_f32(softfloat_propagateNaNF32UI(softfloat_propagateNaNF32UI(uiA, uiB), uiC));
        } else {
            uint32_t const magBits = expA | sigA;
            if (magBits) {
                uint32_t const uiZ = packToF32UI(signProd, 0xFF, 0);
                if (expC != 0xFF) {
                    return ui32_as_f32(uiZ);
                } else if (sigC) {
                    return ui32_as_f32(softfloat_propagateNaNF32UI(uiZ, uiC));
                } else if (signProd == signC) {
                    return ui32_as_f32(uiZ);
                }
            }
            softfloat_raiseFlags(softfloat_flag_invalid);
            return ui32_as_f32(softfloat_propagateNaNF32UI(defaultNaNF32UI, uiC));
        }
    }
    if (expC == 0xFF) {
        return ui32_as_f32(sigC ? softfloat_propagateNaNF32UI(0, uiC) : uiC);
    }

    if (!expA) {
        if (!sigA) {
            return ui32_as_f32(0 == (expC | sigC) && signProd != signC ? packToF32UI(softfloat_roundingMode == softfloat_round_min, 0, 0) : uiC);
        } else {
            struct exp16_sig32 const normExpSig = softfloat_normSubnormalF32Sig(sigA);
            expA = normExpSig.exp;
            sigA = normExpSig.sig;
        }
    }
    if (!expB) {
        if (!sigB) {
            return ui32_as_f32(!(expC | sigC) && signProd != signC ? packToF32UI(softfloat_roundingMode == softfloat_round_min, 0, 0) : uiC);
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
                            return ui32_as_f32(packToF32UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
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

