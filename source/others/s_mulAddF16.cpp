
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

#include "target.hpp"

namespace softfloat {
namespace internals {

float16_t
softfloat_mulAddF16(uint16_t uiA,
                    uint16_t uiB,
                    uint16_t uiC,
                    Mul_add_operations op)
{
    using namespace softfloat::internals;
    bool const signA = signF16UI(uiA);
    int8_t expA = expF16UI(uiA);
    uint16_t sigA = fracF16UI(uiA);
    bool const signB = signF16UI(uiB);
    int8_t expB = expF16UI(uiB);
    uint16_t sigB = fracF16UI(uiB);
    bool signC = signF16UI(uiC) ^ (softfloat_mulAdd_subC == op);
    int8_t expC = expF16UI(uiC);
    uint16_t sigC = fracF16UI(uiC);
    bool signProd = signA ^ signB ^ (softfloat_mulAdd_subProd == op);

    if (expA == 0x1F) {
        if (sigA || ((expB == 0x1F) && sigB)) {
            return u_as_f_16(softfloat_propagateNaNF16UI(softfloat_propagateNaNF16UI(uiA, uiB), uiC));
        }

        uint16_t const magBits = static_cast<uint16_t>(expB | sigB);

        if (magBits) {
            uint16_t const uiZ = packToF16UI(signProd, 0x1F, 0);

            if (expC != 0x1F) {
                return u_as_f_16(uiZ);
            }

            if (sigC) {
                return u_as_f_16(softfloat_propagateNaNF16UI(uiZ, uiC));
            }

            if (signProd == signC) {
                return u_as_f_16(uiZ);
            }
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_16(softfloat_propagateNaNF16UI(defaultNaNF16UI, uiC));
    }

    if (expB == 0x1F) {
        if (sigB) {
            return u_as_f_16(softfloat_propagateNaNF16UI(softfloat_propagateNaNF16UI(uiA, uiB), uiC));
        }

        uint16_t const magBits = static_cast<uint16_t>(expA | sigA);

        if (magBits) {
            uint16_t const uiZ = packToF16UI(signProd, 0x1F, 0);

            if (expC != 0x1F) {
                return u_as_f_16(uiZ);
            }

            if (sigC) {
                return u_as_f_16(softfloat_propagateNaNF16UI(uiZ, uiC));
            }

            if (signProd == signC) {
                return u_as_f_16(uiZ);
            }
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_16(softfloat_propagateNaNF16UI(defaultNaNF16UI, uiC));
    }

    if (expC == 0x1F) {
        return u_as_f_16(sigC ? softfloat_propagateNaNF16UI(0, uiC) : uiC);
    }

    if (!expA) {
        if (!sigA) {
            return u_as_f_16(!(expC | sigC) && (signProd != signC) ? packToF16UI(softfloat_roundingMode == softfloat_round_min, 0, 0) : uiC);
        }

        exp8_sig16 const normExpSig = softfloat_normSubnormalF16Sig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (!expB) {
        if (!sigB) {
            return u_as_f_16(!(expC | sigC) && (signProd != signC) ? packToF16UI(softfloat_roundingMode == softfloat_round_min, 0, 0) : uiC);
        }

        exp8_sig16 const normExpSig = softfloat_normSubnormalF16Sig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int8_t expProd = expA + expB - 0xE;
    sigA = static_cast<uint16_t>((sigA | 0x0400) << 4);
    sigB = static_cast<uint16_t>((sigB | 0x0400) << 4);
    uint32_t sigProd = static_cast<uint32_t>(sigA) * sigB;
    int8_t expZ;
    uint16_t sigZ;
    uint32_t sig32Z;

    if (sigProd < 0x20000000) {
        --expProd;
        sigProd <<= 1;
    }

    bool signZ = signProd;

    if (!expC) {
        if (!sigC) {
            expZ = expProd - 1;
            sigZ = static_cast<uint16_t>(sigProd >> 15 | !!(0 != (sigProd & 0x7FFF)));
            return softfloat_roundPackToF16(signZ, expZ, sigZ);
        }

        exp8_sig16 const normExpSig = softfloat_normSubnormalF16Sig(sigC);
        expC = normExpSig.exp;
        sigC = normExpSig.sig;
    }

    sigC = static_cast<uint16_t>((sigC | 0x0400) << 3);

    int8_t const expDiff = expProd - expC;

    if (signProd == signC) {
        if (expDiff <= 0) {
            expZ = expC;
            /** @todo Warning   C4244   '=': conversion from 'uint32_t' to 'uint16_t', possible loss of data */
            sigZ = static_cast<uint16_t>(sigC + softfloat_shiftRightJam32(sigProd, static_cast<uint16_t>(16 - expDiff)));
        } else {
            expZ = expProd;
            sig32Z =
                sigProd +
                softfloat_shiftRightJam32(static_cast<uint32_t>(sigC) << 16, static_cast<uint16_t>(expDiff));
            sigZ = sig32Z >> 16 | ((sig32Z & 0xFFFF) != 0);
        }

        if (sigZ < 0x4000) {
            --expZ;
            sigZ <<= 1;
        }
    } else {
        uint32_t sig32C = static_cast<uint32_t>(sigC) << 16;

        if (expDiff < 0) {
            signZ = signC;
            expZ = expC;
            sig32Z = sig32C - softfloat_shiftRightJam32(sigProd, static_cast<uint16_t>(-expDiff));
        } else if (!expDiff) {
            expZ = expProd;
            sig32Z = sigProd - sig32C;

            if (!sig32Z) {
                return u_as_f_16(packToF16UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
            }

            if (sig32Z & 0x80000000) {
                signZ = !signZ;
                sig32Z = static_cast<uint32_t>(-static_cast<int32_t>(sig32Z));
            }
        } else {
            expZ = expProd;
            sig32Z = sigProd - softfloat_shiftRightJam32(sig32C, static_cast<uint16_t>(expDiff));
        }

        int8_t shiftDist = softfloat_countLeadingZeros32(sig32Z) - 1;
        expZ -= shiftDist;
        shiftDist -= 16;

        sigZ = static_cast<uint16_t>(
            shiftDist < 0 ?
            sig32Z >> -shiftDist | !!(0 != static_cast<uint32_t>(sig32Z << (shiftDist & 31))) :
            sig32Z << shiftDist);
    }

    return softfloat_roundPackToF16(signZ, expZ, sigZ);
}

}  // namespace internals
}  // namespace softfloat
