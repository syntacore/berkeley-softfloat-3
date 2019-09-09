/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

@copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
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
*/
/*
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

namespace softfloat {
namespace internals {
namespace slow_int64 {

void
softfloat_mulAddF128M(Mul_add_operations op,
                      uint32_t const* const aWPtr,
                      uint32_t const* const bWPtr,
                      uint32_t const* const cWPtr,
                      uint32_t* const zWPtr)
{

    uint32_t const uiA96 = aWPtr[indexWordHi(4)];
    int32_t expA = expF128UI96(uiA96);
    uint32_t const uiB96 = bWPtr[indexWordHi(4)];
    int32_t expB = expF128UI96(uiB96);
    uint32_t const uiC96 = cWPtr[indexWordHi(4)];
    bool const signC = is_sign(uiC96) != (softfloat_mulAdd_subC == op);
    int32_t expC = expF128UI96(uiC96);
    bool signProd = (is_sign(uiA96) != is_sign(uiB96)) != (softfloat_mulAdd_subProd == op);

    bool prodIsInfinite = false;

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (softfloat_tryPropagateNaNF128M(aWPtr, bWPtr, zWPtr)) {
            propagate_NaN_F128M(zWPtr, cWPtr, zWPtr);
            return;
        }

        if (0 == (uiA96 << 1)) {
            if (!(aWPtr[indexWord(4, 2)] | aWPtr[indexWord(4, 1)] | aWPtr[indexWord(4, 0)])) {
                softfloat_invalidF128M(zWPtr);
                propagate_NaN_F128M(zWPtr, cWPtr, zWPtr);
                return;
            }
        } else if (0 == (uiB96 << 1)) {
            if (0 == (bWPtr[indexWord(4, 2)] | bWPtr[indexWord(4, 1)] | bWPtr[indexWord(4, 0)])) {
                softfloat_invalidF128M(zWPtr);
                propagate_NaN_F128M(zWPtr, cWPtr, zWPtr);
                return;
            }
        }

        prodIsInfinite = true;
    }

    if (0x7FFF == expC) {
        if (fracF128UI96(uiC96) || (cWPtr[indexWord(4, 2)] | cWPtr[indexWord(4, 1)] | cWPtr[indexWord(4, 0)])) {
            zWPtr[indexWordHi(4)] = 0;
            propagate_NaN_F128M(zWPtr, cWPtr, zWPtr);
            return;
        }

        if (prodIsInfinite && (signProd != signC)) {
            softfloat_invalidF128M(zWPtr);
            propagate_NaN_F128M(zWPtr, cWPtr, zWPtr);
            return;
        }

        zWPtr[indexWordHi(4)] = uiC96;
        zWPtr[indexWord(4, 2)] = cWPtr[indexWord(4, 2)];
        zWPtr[indexWord(4, 1)] = cWPtr[indexWord(4, 1)];
        zWPtr[indexWord(4, 0)] = cWPtr[indexWord(4, 0)];
        return;
    }

    if (prodIsInfinite) {
        zWPtr[indexWordHi(4)] = packToF128UI96(signProd, 0x7FFF, 0);
        zWPtr[indexWord(4, 2)] = 0;
        zWPtr[indexWord(4, 1)] = 0;
        zWPtr[indexWord(4, 0)] = 0;
        return;
    }

    uint32_t sigA[4];

    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    if (0 == expA) {
        sigA[indexWordHi(4)] = fracF128UI96(uiA96) | 0x00010000;
        sigA[indexWord(4, 2)] = aWPtr[indexWord(4, 2)];
        sigA[indexWord(4, 1)] = aWPtr[indexWord(4, 1)];
        sigA[indexWord(4, 0)] = aWPtr[indexWord(4, 0)];
    } else {
        expA = softfloat_shiftNormSigF128M(aWPtr, 0, sigA);

        if (expA == -128) {
            if (!static_cast<uint32_t>(uiC96 << 1) && (signProd != signC) && !cWPtr[indexWord(4, 2)] && !(cWPtr[indexWord(4, 1)] | cWPtr[indexWord(4, 0)])) {
                zWPtr[indexWordHi(4)] = packToF128UI96(softfloat_round_min == softfloat_roundingMode, 0u, 0u);
                zWPtr[indexWord(4, 2)] = 0;
                zWPtr[indexWord(4, 1)] = 0;
                zWPtr[indexWord(4, 0)] = 0;
                return;
            }

            zWPtr[indexWordHi(4)] = uiC96;
            zWPtr[indexWord(4, 2)] = cWPtr[indexWord(4, 2)];
            zWPtr[indexWord(4, 1)] = cWPtr[indexWord(4, 1)];
            zWPtr[indexWord(4, 0)] = cWPtr[indexWord(4, 0)];
            return;
        }
    }

    uint32_t sigX[5];

    if (expB) {
        sigX[indexWordHi(4)] = fracF128UI96(uiB96) | 0x00010000;
        sigX[indexWord(4, 2)] = bWPtr[indexWord(4, 2)];
        sigX[indexWord(4, 1)] = bWPtr[indexWord(4, 1)];
        sigX[indexWord(4, 0)] = bWPtr[indexWord(4, 0)];
    } else {
        expB = softfloat_shiftNormSigF128M(bWPtr, 0, sigX);

        if (expB == -128) {
            if (!static_cast<uint32_t>(uiC96 << 1) && (signProd != signC) && !cWPtr[indexWord(4, 2)] && !(cWPtr[indexWord(4, 1)] | cWPtr[indexWord(4, 0)])) {
                zWPtr[indexWordHi(4)] = packToF128UI96(softfloat_round_min == softfloat_roundingMode, 0u, 0u);
                zWPtr[indexWord(4, 2)] = 0;
                zWPtr[indexWord(4, 1)] = 0;
                zWPtr[indexWord(4, 0)] = 0;
                return;
            }

            zWPtr[indexWordHi(4)] = uiC96;
            zWPtr[indexWord(4, 2)] = cWPtr[indexWord(4, 2)];
            zWPtr[indexWord(4, 1)] = cWPtr[indexWord(4, 1)];
            zWPtr[indexWord(4, 0)] = cWPtr[indexWord(4, 0)];
            return;
        }
    }

    int32_t const expProd = expA + expB - 0x3FF0;
    uint32_t sigProd[8];
    softfloat_mul128MTo256M(sigA, sigX, sigProd);

    uint32_t wordSig = fracF128UI96(uiC96);

    if (expC) {
        --expC;
        wordSig |= 0x00010000;
    }

    sigX[indexWordHi(5)] = wordSig;
    sigX[indexWord(5, 3)] = cWPtr[indexWord(4, 2)];
    sigX[indexWord(5, 2)] = cWPtr[indexWord(4, 1)];
    sigX[indexWord(5, 1)] = cWPtr[indexWord(4, 0)];

    bool doSub = signProd != signC;
    auto const addCarryMRoutinePtr = doSub ? softfloat_addComplCarryM : softfloat_addCarryM;
    int32_t expDiff = expProd - expC;
    uint32_t* extSigPtr;

    bool signZ;
    bool carry;
    int32_t expZ;

    if (expDiff <= 0) {
        signZ = signC;
        expZ = expC;

        if (sigProd[indexWord(8, 2)] || (sigProd[indexWord(8, 1)] | sigProd[indexWord(8, 0)])) {
            sigProd[indexWord(8, 3)] |= 1;
        }

        extSigPtr = &sigProd[indexMultiwordHi(8, 5)];

        if (expDiff) {
            softfloat_shiftRightJam160M(extSigPtr, static_cast<uint8_t>(-expDiff), extSigPtr);
        }

        carry = false;

        if (doSub) {
            wordSig = extSigPtr[indexWordLo(5)];
            extSigPtr[indexWordLo(5)] = static_cast<uint32_t>(-static_cast<int32_t>(wordSig));
            carry = 0 == wordSig;
        }

        (*addCarryMRoutinePtr)(
            4u,
            &sigX[indexMultiwordHi(5, 4)],
            extSigPtr + indexMultiwordHi(5, 4),
            carry,
            extSigPtr + indexMultiwordHi(5, 4)
        );
        wordSig = extSigPtr[indexWordHi(5)];

        if (!expZ) {
            if (wordSig & 0x80000000) {
                signZ = !signZ;
                softfloat_negX160M(extSigPtr);
                wordSig = extSigPtr[indexWordHi(5)];
            }

            if (wordSig || (extSigPtr[indexWord(5, 3)] | extSigPtr[indexWord(5, 2)]) || (extSigPtr[indexWord(5, 1)] | extSigPtr[indexWord(5, 0)])) {
                if (wordSig < 0x00010000) {
                    softfloat_normRoundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
                    return;
                }

                if (0x00020000 <= wordSig) {
                    ++expZ;
                    softfloat_shortShiftRightJam160M(extSigPtr, 1, extSigPtr);
                }

                softfloat_roundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
                return;
            }

            zWPtr[indexWordHi(4)] = packToF128UI96(softfloat_round_min == softfloat_roundingMode, 0u, 0u);
            zWPtr[indexWord(4, 2)] = 0;
            zWPtr[indexWord(4, 1)] = 0;
            zWPtr[indexWord(4, 0)] = 0;
            return;
        }

        if (wordSig < 0x00010000) {
            --expZ;
            softfloat_add160M(extSigPtr, extSigPtr, extSigPtr);
            softfloat_roundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
            return;
        }

        if (0x00020000 <= wordSig) {
            ++expZ;
            softfloat_shortShiftRightJam160M(extSigPtr, 1, extSigPtr);
        }

        softfloat_roundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
        return;
    }

    signZ = signProd;
    expZ = expProd;
    sigX[indexWordLo(5)] = 0;
    expDiff -= 128;

    uint32_t* ptr;

    if (0 <= expDiff) {
        if (expDiff) {
            softfloat_shiftRightJam160M(sigX, static_cast<uint8_t>(expDiff), sigX);
        }

        wordSig = sigX[indexWordLo(5)];
        carry = 0;

        if (doSub) {
            carry = 0 == wordSig;
            wordSig = static_cast<uint32_t>(-static_cast<int32_t>(wordSig));
        }

        carry =
            (*addCarryMRoutinePtr)(
                4,
                &sigProd[indexMultiwordLo(8, 4)],
                &sigX[indexMultiwordHi(5, 4)],
                carry,
                &sigProd[indexMultiwordLo(8, 4)]
            );
        sigProd[indexWord(8, 2)] |= wordSig;
        ptr = &sigProd[indexWord(8, 4)];
    } else {
        uint8_t const shiftDist = static_cast<uint8_t>(expDiff & 31);

        if (shiftDist) {
            softfloat_shortShiftRight160M(sigX, shiftDist, sigX);
        }

        //TODO: V610 Unspecified behavior. Check the shift operator '>>='. The left operand is negative ('expDiff' = [-127..-1]). https://www.viva64.com/en/w/v610/
        expDiff >>= 5;
        extSigPtr = &sigProd[indexMultiwordLo(8, 5)] - wordIncr + expDiff * -wordIncr;
        carry = (*addCarryMRoutinePtr)(5, extSigPtr, sigX, doSub, extSigPtr);

        if (expDiff == -4) {
            wordSig = sigProd[indexWordHi(8)];

            if (wordSig & 0x80000000) {
                signZ = !signZ;
                softfloat_negX256M(sigProd);
                wordSig = sigProd[indexWordHi(8)];
            }

            if (wordSig) {
                if (sigProd[indexWord(8, 2)] || (sigProd[indexWord(8, 1)] | sigProd[indexWord(8, 0)])) {
                    sigProd[indexWord(8, 3)] |= 1;
                }

                extSigPtr = &sigProd[indexMultiwordHi(8, 5)];
                wordSig = extSigPtr[indexWordHi(5)];

                if (wordSig < 0x00010000) {
                    softfloat_normRoundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
                    return;
                }

                if (0x00020000 <= wordSig) {
                    ++expZ;
                    softfloat_shortShiftRightJam160M(extSigPtr, 1, extSigPtr);
                }

                softfloat_roundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
                return;
            }

            wordSig = sigProd[indexWord(8, 6)];

            if (0x00040000 <= wordSig) {
                if (sigProd[indexWord(8, 2)] || (sigProd[indexWord(8, 1)] | sigProd[indexWord(8, 0)])) {
                    sigProd[indexWord(8, 3)] |= 1;
                }

                extSigPtr = &sigProd[indexMultiwordHi(8, 5)];
                wordSig = extSigPtr[indexWordHi(5)];

                if (wordSig < 0x00010000) {
                    softfloat_normRoundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
                    return;
                }

                if (0x00020000 <= wordSig) {
                    ++expZ;
                    softfloat_shortShiftRightJam160M(extSigPtr, 1, extSigPtr);
                }

                softfloat_roundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
                return;
            }

            expZ -= 32;
            extSigPtr = &sigProd[indexMultiwordHi(8, 5)] - wordIncr;

            for (;;) {
                if (wordSig) {
                    break;
                }

                wordSig = extSigPtr[indexWord(5, 3)];

                if (0x00040000 <= wordSig) {
                    break;
                }

                expZ -= 32;
                extSigPtr -= wordIncr;

                if (extSigPtr == &sigProd[indexMultiwordLo(8, 5)]) {
                    if (wordSig || (extSigPtr[indexWord(5, 3)] | extSigPtr[indexWord(5, 2)]) || (extSigPtr[indexWord(5, 1)] | extSigPtr[indexWord(5, 0)])) {
                        if (wordSig < 0x00010000) {
                            softfloat_normRoundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
                            return;
                        }

                        if (0x00020000 <= wordSig) {
                            ++expZ;
                            softfloat_shortShiftRightJam160M(extSigPtr, 1, extSigPtr);
                        }

                        softfloat_roundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
                        return;
                    }

                    zWPtr[indexWordHi(4)] = packToF128UI96(softfloat_round_min == softfloat_roundingMode, 0u, 0u);
                    zWPtr[indexWord(4, 2)] = 0;
                    zWPtr[indexWord(4, 1)] = 0;
                    zWPtr[indexWord(4, 0)] = 0;
                    return;
                }
            }

            ptr = extSigPtr + indexWordLo(5);

            do {
                ptr -= wordIncr;

                if (*ptr) {
                    extSigPtr[indexWordLo(5)] |= 1;
                    break;
                }
            } while (ptr != &sigProd[indexWordLo(8)]);

            wordSig = extSigPtr[indexWordHi(5)];

            if (wordSig < 0x00010000) {
                softfloat_normRoundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
                return;
            }

            if (0x00020000 <= wordSig) {
                ++expZ;
                softfloat_shortShiftRightJam160M(extSigPtr, 1, extSigPtr);
            }

            softfloat_roundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
            return;
        }

        ptr = extSigPtr + indexWordHi(5) + wordIncr;
    }

    if (carry != doSub) {
        if (doSub) {
            do {
                wordSig = *ptr;
                *ptr = wordSig - 1;
                ptr += wordIncr;
            } while (!wordSig);
        } else {
            do {
                wordSig = *ptr + 1;
                *ptr = wordSig;
                ptr += wordIncr;
            } while (!wordSig);
        }
    }

    if (sigProd[indexWord(8, 2)] || (sigProd[indexWord(8, 1)] | sigProd[indexWord(8, 0)])) {
        sigProd[indexWord(8, 3)] |= 1;
    }

    extSigPtr = &sigProd[indexMultiwordHi(8, 5)];
    wordSig = extSigPtr[indexWordHi(5)];

    if (wordSig < 0x00010000) {
        softfloat_normRoundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
        return;
    }

    if (0x00020000 <= wordSig) {
        ++expZ;
        softfloat_shortShiftRightJam160M(extSigPtr, 1, extSigPtr);
    }

    softfloat_roundPackMToF128M(signZ, expZ, extSigPtr, zWPtr);
}

}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
