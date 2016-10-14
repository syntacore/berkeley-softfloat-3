
/*============================================================================

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

=============================================================================*/

#include "internals.h"

#include "specialize.h"
#include "softfloat/functions.h"

float16_t
softfloat_mulAddF16(uint16_t uiA, uint16_t uiB, uint16_t uiC, uint8_t op)
{
    uint16_t magBits, uiZ;
    struct exp8_sig16 normExpSig;
    int8_t expProd;
    uint32_t sigProd;
    bool signZ;
    int8_t expZ;
    uint16_t sigZ;
    int8_t expDiff;
    uint32_t sig32Z, sig32C;
    int8_t shiftDist;
    union ui16_f16 uZ;


    bool const signA = signF16UI(uiA);
    int8_t expA = expF16UI(uiA);
    uint16_t sigA = fracF16UI(uiA);
    bool const signB = signF16UI(uiB);
    int8_t expB = expF16UI(uiB);
    uint16_t sigB = fracF16UI(uiB);
    bool signC = signF16UI(uiC) ^ (op == softfloat_mulAdd_subC);
    int8_t expC = expF16UI(uiC);
    uint16_t sigC = fracF16UI(uiC);
    bool signProd = signA ^ signB ^ (op == softfloat_mulAdd_subProd);

    if (expA == 0x1F) {
        if (sigA || ((expB == 0x1F) && sigB)) {
            goto propagateNaN_ABC;
        } else {
            magBits = expB | sigB;
            goto infProdArg;
        }
    }
    if (expB == 0x1F) {
        if (sigB) {
            goto propagateNaN_ABC;
        } else {
            magBits = expA | sigA;
            goto infProdArg;
        }
    }
    if (expC == 0x1F) {
        if (sigC) {
            uiZ = 0;
            goto propagateNaN_ZC;
        } else {
            uiZ = uiC;
            goto uiZ;
        }
    }
    if (!expA) {
        if (!sigA) goto zeroProd;
        normExpSig = softfloat_normSubnormalF16Sig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }
    if (!expB) {
        if (!sigB) goto zeroProd;
        normExpSig = softfloat_normSubnormalF16Sig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    expProd = expA + expB - 0xE;
    sigA = (sigA | 0x0400) << 4;
    sigB = (sigB | 0x0400) << 4;
    sigProd = (uint32_t)sigA * sigB;
    if (sigProd < 0x20000000) {
        --expProd;
        sigProd <<= 1;
    }
    signZ = signProd;
    if (!expC) {
        if (!sigC) {
            expZ = expProd - 1;
            sigZ = sigProd >> 15 | ((sigProd & 0x7FFF) != 0);
            goto roundPack;
        }
        normExpSig = softfloat_normSubnormalF16Sig(sigC);
        expC = normExpSig.exp;
        sigC = normExpSig.sig;
    }
    sigC = (sigC | 0x0400) << 3;
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
    expDiff = expProd - expC;
    if (signProd == signC) {
        /*--------------------------------------------------------------------
        *--------------------------------------------------------------------*/
        if (expDiff <= 0) {
            expZ = expC;
            sigZ = sigC + softfloat_shiftRightJam32(sigProd, 16 - expDiff);
        } else {
            expZ = expProd;
            sig32Z =
                sigProd
                + softfloat_shiftRightJam32(
                (uint32_t)sigC << 16, expDiff);
            sigZ = sig32Z >> 16 | ((sig32Z & 0xFFFF) != 0);
        }
        if (sigZ < 0x4000) {
            --expZ;
            sigZ <<= 1;
        }
    } else {
        /*--------------------------------------------------------------------
        *--------------------------------------------------------------------*/
        sig32C = (uint32_t)sigC << 16;
        if (expDiff < 0) {
            signZ = signC;
            expZ = expC;
            sig32Z = sig32C - softfloat_shiftRightJam32(sigProd, -expDiff);
        } else if (!expDiff) {
            expZ = expProd;
            sig32Z = sigProd - sig32C;
            if (!sig32Z) goto completeCancellation;
            if (sig32Z & 0x80000000) {
                signZ = !signZ;
                sig32Z = -(int32_t)sig32Z;
            }
        } else {
            expZ = expProd;
            sig32Z = sigProd - softfloat_shiftRightJam32(sig32C, expDiff);
        }
        shiftDist = softfloat_countLeadingZeros32(sig32Z) - 1;
        expZ -= shiftDist;
        shiftDist -= 16;
        if (shiftDist < 0) {
            sigZ =
                sig32Z >> (-shiftDist)
                | ((uint32_t)(sig32Z << (shiftDist & 31)) != 0);
        } else {
            sigZ = (uint16_t)sig32Z << shiftDist;
        }
    }
roundPack:
    return softfloat_roundPackToF16(signZ, expZ, sigZ);

propagateNaN_ABC:
    uiZ = softfloat_propagateNaNF16UI(uiA, uiB);
    goto propagateNaN_ZC;
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
infProdArg:
    if (magBits) {
        uiZ = packToF16UI(signProd, 0x1F, 0);
        if (expC != 0x1F) {
            goto uiZ;
        } else if (sigC) {
            goto propagateNaN_ZC;
        } else if (signProd == signC) {
            goto uiZ;
        }
    }
    softfloat_raiseFlags(softfloat_flag_invalid);
    uiZ = defaultNaNF16UI;
propagateNaN_ZC:
    uiZ = softfloat_propagateNaNF16UI(uiZ, uiC);
    goto uiZ;

zeroProd:
    uiZ = uiC;
    if (!(expC | sigC) && (signProd != signC)) {
completeCancellation:
        uiZ = packToF16UI(softfloat_roundingMode == softfloat_round_min, 0, 0);
    }
uiZ:
    uZ.ui = uiZ;
    return uZ.f;
}
