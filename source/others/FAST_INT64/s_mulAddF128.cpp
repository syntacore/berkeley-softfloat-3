
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

namespace softfloat {
namespace internals {

float128_t
softfloat_mulAddF128(uint64_t const uiA64,
                     uint64_t const uiA0,
                     uint64_t const uiB64,
                     uint64_t const uiB0,
                     uint64_t const uiC64,
                     uint64_t const uiC0,
                     Mul_add_operations const op)
{
    bool const signA = is_sign(uiA64);
    int32_t expA = expF128UI64(uiA64);
    uint128 sigA{fracF128UI64(uiA64), uiA0};
    bool const signB = is_sign(uiB64);
    int32_t expB = expF128UI64(uiB64);
    uint128 sigB{fracF128UI64(uiB64), uiB0};
    bool const signC = is_sign(uiC64) != (softfloat_mulAdd_subC == op);
    int32_t expC = expF128UI64(uiC64);
    uint128 sigC{fracF128UI64(uiC64), uiC0};
    bool signZ = (signA != signB) != (softfloat_mulAdd_subProd == op);

    if (0x7FFF == expA) {
        if (0 != (sigA.v64 | sigA.v0) || (expB == 0x7FFF && 0 != (sigB.v64 | sigB.v0))) {
            auto const uiZ_1 = softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0);
            return static_cast<float128_t>(softfloat_propagateNaNF128UI(uiZ_1.v64, uiZ_1.v0, uiC64, uiC0));
        }

        if (0 != (expB | sigB.v64 | sigB.v0)) {
            uint128 const uiZ_1 = uint128{packToF128UI64(signZ, 0x7FFF, 0), 0};

            if (0x7FFF != expC) {
                return static_cast<float128_t>(uiZ_1);
            }

            if (0 != (sigC.v64 | sigC.v0)) {
                return static_cast<float128_t>(softfloat_propagateNaNF128UI(uiZ_1.v64, uiZ_1.v0, uiC64, uiC0));
            }

            if (signZ == signC) {
                return static_cast<float128_t>(uiZ_1);
            }
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return static_cast<float128_t>(softfloat_propagateNaNF128UI(defaultNaNF128UI64, defaultNaNF128UI0, uiC64, uiC0));
    }

    if (0x7FFF == expB) {
        if (0 != (sigB.v64 | sigB.v0)) {
            auto const uiZ_1 = softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0);
            return static_cast<float128_t>(softfloat_propagateNaNF128UI(uiZ_1.v64, uiZ_1.v0, uiC64, uiC0));
        }

        if (0 != (expA | sigA.v64 | sigA.v0)) {
            auto const uiZ_1 = uint128{packToF128UI64(signZ, 0x7FFF, 0), 0};

            if (0x7FFF != expC) {
                return static_cast<float128_t>(uiZ_1);
            }

            if (0 != (sigC.v64 | sigC.v0)) {
                return static_cast<float128_t>(softfloat_propagateNaNF128UI(uiZ_1.v64, uiZ_1.v0, uiC64, uiC0));
            }

            if (signZ == signC) {
                return static_cast<float128_t>(uiZ_1);
            }
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return static_cast<float128_t>(softfloat_propagateNaNF128UI(defaultNaNF128UI64, defaultNaNF128UI0, uiC64, uiC0));
    }

    if (0x7FFF == expC) {
        if (0 != (sigC.v64 | sigC.v0)) {
            return static_cast<float128_t>(softfloat_propagateNaNF128UI(0, 0, uiC64, uiC0));
        }

        return static_cast<float128_t>(uint128{uiC64, uiC0});
    }

    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    if (0 == expA) {
        if (0 == (sigA.v64 | sigA.v0)) {
            if (0 == (expC | sigC.v64 | sigC.v0) && signZ != signC) {
                return static_cast<float128_t>(uint128{packToF128UI64((softfloat_roundingMode == softfloat_round_min), 0, 0), 0});
            }

            return static_cast<float128_t>(uint128{uiC64, uiC0});
        }

        exp32_sig128 const normExpSig = softfloat_normSubnormalF128Sig(sigA.v64, sigA.v0);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (0 == expB) {
        if (0 == (sigB.v64 | sigB.v0)) {
            if (0 == (expC | sigC.v64 | sigC.v0) && signZ != signC) {
                return static_cast<float128_t>(uint128{packToF128UI64((softfloat_roundingMode == softfloat_round_min), 0, 0), 0});
            }

            return static_cast<float128_t>(uint128{uiC64, uiC0});
        }

        exp32_sig128 const normExpSig = softfloat_normSubnormalF128Sig(sigB.v64, sigB.v0);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int32_t expZ = expA + expB - 0x3FFE;
    sigA.v64 |= UINT64_C(0x0001000000000000);
    sigB.v64 |= UINT64_C(0x0001000000000000);
    sigA = softfloat_shortShiftLeft128(sigA.v64, sigA.v0, 8);
    sigB = softfloat_shortShiftLeft128(sigB.v64, sigB.v0, 15);

    uint64_t sig256Z[4];
    softfloat_mul128To256M(sigA.v64, sigA.v0, sigB.v64, sigB.v0, sig256Z);

    uint128 sigZ;
    sigZ.v64 = sig256Z[indexWord(4, 3)];
    sigZ.v0 = sig256Z[indexWord(4, 2)];

    int32_t shiftDist_1 = 0;

    if (0 == (sigZ.v64 & UINT64_C(0x0100000000000000))) {
        --expZ;
        shiftDist_1 = -1;
    }

    if (0 == expC) {
        if (0 == (sigC.v64 | sigC.v0)) {
            shiftDist_1 += 8;
            auto const sigZExtra_2 = sig256Z[indexWord(4, 1)] | sig256Z[indexWord(4, 0)];
            auto const sigZExtra_1 = static_cast<uint64_t>(sigZ.v0 << (64 - shiftDist_1)) | (sigZExtra_2 != 0);
            auto const sigZ_1 = softfloat_shortShiftRight128(sigZ.v64, sigZ.v0, static_cast<uint8_t>(shiftDist_1));
            return softfloat_roundPackToF128(signZ, expZ - 1, sigZ_1.v64, sigZ_1.v0, sigZExtra_1);
        }

        exp32_sig128 const normExpSig = softfloat_normSubnormalF128Sig(sigC.v64, sigC.v0);
        expC = normExpSig.exp;
        sigC = normExpSig.sig;
    }

    sigC.v64 |= UINT64_C(0x0001000000000000);
    sigC = softfloat_shortShiftLeft128(sigC.v64, sigC.v0, 8);

    int32_t const expDiff = expZ - expC;
    uint64_t sig256C[4];

    if (expDiff < 0) {
        expZ = expC;

        if (signZ == signC || expDiff < -1) {
            shiftDist_1 -= expDiff;

            if (0 != shiftDist_1) {
                sigZ = softfloat_shiftRightJam128(sigZ.v64, sigZ.v0, static_cast<uint32_t>(shiftDist_1));
            }
        } else {
            if (0 == shiftDist_1) {
                uint128 const x128 = softfloat_shortShiftRight128(sig256Z[indexWord(4, 1)], sig256Z[indexWord(4, 0)], 1);
                sig256Z[indexWord(4, 1)] = (sigZ.v0 << 63) | x128.v64;
                sig256Z[indexWord(4, 0)] = x128.v0;
                sigZ = softfloat_shortShiftRight128(sigZ.v64, sigZ.v0, 1);
                sig256Z[indexWord(4, 3)] = sigZ.v64;
                sig256Z[indexWord(4, 2)] = sigZ.v0;
            }
        }
    } else {
        if (0 != shiftDist_1) {
            softfloat_add256M(sig256Z, sig256Z, sig256Z);
        }

        if (0 == expDiff) {
            sigZ.v64 = sig256Z[indexWord(4, 3)];
            sigZ.v0 = sig256Z[indexWord(4, 2)];
        } else {
            sig256C[indexWord(4, 3)] = sigC.v64;
            sig256C[indexWord(4, 2)] = sigC.v0;
            sig256C[indexWord(4, 1)] = 0;
            sig256C[indexWord(4, 0)] = 0;
            softfloat_shiftRightJam256M(sig256C, static_cast<uint32_t>(expDiff), sig256C);
        }
    }

    {
        int32_t shiftDist_2 = 8;

        if (signZ == signC) {
            if (expDiff <= 0) {
                sigZ = softfloat_add128(sigC.v64, sigC.v0, sigZ.v64, sigZ.v0);
            } else {
                softfloat_add256M(sig256Z, sig256C, sig256Z);
                sigZ.v64 = sig256Z[indexWord(4, 3)];
                sigZ.v0 = sig256Z[indexWord(4, 2)];
            }

            if (sigZ.v64 & UINT64_C(0x0200000000000000)) {
                ++expZ;
                shiftDist_2 = 9;
            }

            auto const sigZExtra_2 = sig256Z[indexWord(4, 1)] | sig256Z[indexWord(4, 0)];
            auto const sigZExtra_1 = static_cast<uint64_t>(sigZ.v0 << (64 - shiftDist_2)) | (sigZExtra_2 != 0);
            auto const sigZ_1 = softfloat_shortShiftRight128(sigZ.v64, sigZ.v0, static_cast<uint8_t>(shiftDist_2));
            return softfloat_roundPackToF128(signZ, expZ - 1, sigZ_1.v64, sigZ_1.v0, sigZExtra_1);
        }

        {
            uint64_t sigZExtra;

            if (expDiff < 0) {
                signZ = signC;

                if (expDiff < -1) {
                    sigZ = softfloat_sub128(sigC.v64, sigC.v0, sigZ.v64, sigZ.v0);
                    sigZExtra = sig256Z[indexWord(4, 1)] | sig256Z[indexWord(4, 0)];

                    if (0 != sigZExtra) {
                        sigZ = softfloat_sub128(sigZ.v64, sigZ.v0, 0, 1);
                    }

                    if (0 == (sigZ.v64 & UINT64_C(0x0100000000000000))) {
                        --expZ;
                        shiftDist_2 = 7;
                    }

                    auto const sigZExtra_1 = static_cast<uint64_t>(sigZ.v0 << (64 - shiftDist_2)) | !!(0 != sigZExtra);
                    auto const sigZ_1 = softfloat_shortShiftRight128(sigZ.v64, sigZ.v0, static_cast<uint8_t>(shiftDist_2));
                    return softfloat_roundPackToF128(signZ, expZ - 1, sigZ_1.v64, sigZ_1.v0, sigZExtra_1);
                }

                sig256C[indexWord(4, 3)] = sigC.v64;
                sig256C[indexWord(4, 2)] = sigC.v0;
                sig256C[indexWord(4, 1)] = 0;
                sig256C[indexWord(4, 0)] = 0;
                softfloat_sub256M(sig256C, sig256Z, sig256Z);
            } else if (0 == expDiff) {
                sigZ = softfloat_sub128(sigZ.v64, sigZ.v0, sigC.v64, sigC.v0);

                if (0 == (sigZ.v64 | sigZ.v0) && 0 == sig256Z[indexWord(4, 1)] && 0 == sig256Z[indexWord(4, 0)]) {
                    return static_cast<float128_t>(uint128{packToF128UI64(softfloat_round_min == softfloat_roundingMode, 0, 0), 0});
                }

                sig256Z[indexWord(4, 3)] = sigZ.v64;
                sig256Z[indexWord(4, 2)] = sigZ.v0;

                if (0 != (sigZ.v64 & UINT64_C(0x8000000000000000))) {
                    signZ = !signZ;
                    uint64_t zero256[4] = INIT_UINTM4(0, 0, 0, 0);
                    softfloat_sub256M(zero256, sig256Z, sig256Z);
                }
            } else {
                softfloat_sub256M(sig256Z, sig256C, sig256Z);

                if (1 < expDiff) {
                    sigZ.v64 = sig256Z[indexWord(4, 3)];
                    sigZ.v0 = sig256Z[indexWord(4, 2)];

                    if (0 == (sigZ.v64 & UINT64_C(0x0100000000000000))) {
                        --expZ;
                        shiftDist_2 = 7;
                    }

                    auto const sigZExtra_2 = sig256Z[indexWord(4, 1)] | sig256Z[indexWord(4, 0)];
                    auto const sigZExtra_1 = static_cast<uint64_t>(sigZ.v0 << (64 - shiftDist_2)) | (sigZExtra_2 != 0);
                    auto const sigZ_1 = softfloat_shortShiftRight128(sigZ.v64, sigZ.v0, static_cast<uint8_t>(shiftDist_2));
                    return softfloat_roundPackToF128(signZ, expZ - 1, sigZ_1.v64, sigZ_1.v0, sigZExtra_1);
                }
            }
        }

        sigZ.v64 = sig256Z[indexWord(4, 3)];
        sigZ.v0 = sig256Z[indexWord(4, 2)];

        {
            uint64_t sigZExtra = sig256Z[indexWord(4, 1)];
            uint64_t const sig256Z0 = sig256Z[indexWord(4, 0)];

            if (0 != sigZ.v64) {
                if (0 != sig256Z0) {
                    sigZExtra |= 1;
                }
            } else {
                expZ -= 64;
                sigZ.v64 = sigZ.v0;
                sigZ.v0 = sigZExtra;
                sigZExtra = sig256Z0;

                if (0 == sigZ.v64) {
                    expZ -= 64;
                    sigZ.v64 = sigZ.v0;
                    sigZ.v0 = sigZExtra;
                    sigZExtra = 0;

                    if (0 == sigZ.v64) {
                        expZ -= 64;
                        sigZ.v64 = sigZ.v0;
                        sigZ.v0 = 0;
                    }
                }
            }

            shiftDist_2 = count_leading_zeros(sigZ.v64);
            expZ += 7 - shiftDist_2;
            shiftDist_2 = 15 - shiftDist_2;

            if (0 < shiftDist_2) {
                auto const sigZExtra_1 = static_cast<uint64_t>(sigZ.v0 << (64 - shiftDist_2)) | !!(0 != sigZExtra);
                auto const sigZ_1 = softfloat_shortShiftRight128(sigZ.v64, sigZ.v0, static_cast<uint8_t>(shiftDist_2));
                return softfloat_roundPackToF128(signZ, expZ - 1, sigZ_1.v64, sigZ_1.v0, sigZExtra_1);
            }

            if (0 != shiftDist_2) {
                shiftDist_2 = -shiftDist_2;
                sigZ = softfloat_shortShiftLeft128(sigZ.v64, sigZ.v0, static_cast<uint8_t>(shiftDist_2));
                uint128 const x128 = softfloat_shortShiftLeft128(0, sigZExtra, static_cast<uint8_t>(shiftDist_2));
                sigZ.v0 |= x128.v64;
                sigZExtra = x128.v0;
            }

            return softfloat_roundPackToF128(signZ, expZ - 1, sigZ.v64, sigZ.v0, sigZExtra);
        }
    }
}

}  // namespace internals
}  // namespace softfloat
