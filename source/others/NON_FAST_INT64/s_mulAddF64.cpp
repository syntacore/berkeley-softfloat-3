
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

#ifdef SOFTFLOAT_FAST_INT64
#error For non-fast int64_t only
#endif

namespace softfloat {
namespace internals {

namespace {

static inline float64_t
mulAdd(Mul_add_operations op,
       uint64_t const uiA,
       uint64_t const uiB,
       uint64_t const uiC
      )
{
    if (softfloat_isSigNaNF64UI(uiC) || softfloat_isNaNF64UI(uiA) || softfloat_isNaNF64UI(uiB)) {
        return u_as_f_64(softfloat_propagateNaNF64UI(softfloat_propagateNaNF64UI(uiA, uiB), uiC));
    }

    bool const signA = is_sign(uiA);
    uint64_t sigA = fracF64UI(uiA);
    bool const signB = is_sign(uiB);
    uint64_t sigB = fracF64UI(uiB);
    bool const signC = is_sign(uiC) != (softfloat_mulAdd_subC == op);
    int16_t expC = expF64UI(uiC);
    uint64_t sigC = fracF64UI(uiC);
    bool signZ = (signA != signB) != (softfloat_mulAdd_subProd == op);

    if (isInf64UI(uiA) || isInf64UI(uiB)) {
        /* a or b is inf, product is inf or undefined, check other operand for zero */
        bool const is_product_undefined = isInf64UI(uiA) ? isZero64UI(uiB) : isZero64UI(uiA);

        if (is_product_undefined) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return u_as_f_64(softfloat_propagateNaNF64UI(defaultNaNF64UI, uiC));
        }

        /* product is inf */
        uint64_t const uiZ = packToF64UI(signZ, 0x7FF, 0);

        if (expC != 0x7FF) {
            /* summand c is finite, return product as result */
            return u_as_f_64(uiZ);
        }

        /* if summand is inf, check for same sign */
        if (!softfloat_isNaNF64UI(uiC) && signZ == signC) {
            /* summands are same sign inf */
            return u_as_f_64(uiZ);
        }

        /* summands are different sign inf or NaN, undefined sum */
        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f_64(softfloat_propagateNaNF64UI(defaultNaNF64UI, uiC));
    }

    if (softfloat_isNaNF64UI(uiC)) {
        return u_as_f_64(softfloat_propagateNaNF64UI(defaultNaNF64UI, uiC));
    }

    if (expC == 0x7FF) {
        /** if c is infinity while a and b are finite, return c */
        return u_as_f_64(uiC);
    }

    int16_t expA = expF64UI(uiA);

    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    if (0 == expA) {
        if (!sigA) {
            return u_as_f_64(0 == (expC | sigC) && signZ != signC ? packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0) : uiC);
        }

        exp16_sig64 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int16_t expB = expF64UI(uiB);

    if (!expB) {
        if (!sigB) {
            return u_as_f_64(0 == (expC | sigC) && signZ != signC ? packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0) : uiC);
        }

        exp16_sig64 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int16_t expZ = expA + expB - 0x3FE;
    sigA = (sigA | UINT64_C(0x0010000000000000)) << 10;
    sigB = (sigB | UINT64_C(0x0010000000000000)) << 11;
    uint32_t sig128Z[4];
    softfloat_mul64To128M(sigA, sigB, sig128Z);
    uint64_t sigZ = static_cast<uint64_t>(sig128Z[indexWord(4, 3)]) << 32 | sig128Z[indexWord(4, 2)];
    int16_t shiftDist = 0;

    if (0 == (sigZ & UINT64_C(0x4000000000000000))) {
        --expZ;
        shiftDist = -1;
    }

    if (0 == expC) {
        if (0 == sigC) {
            if (0 != shiftDist) {
                sigZ <<= 1;
            }

            return
                softfloat_roundPackToF64(signZ, expZ - 1,
                                         sigZ | (sig128Z[indexWord(4, 1)] || sig128Z[indexWord(4, 0)]));
        }

        exp16_sig64 const normExpSig(sigC);
        expC = normExpSig.exp;
        sigC = normExpSig.sig;
    }

    sigC = (sigC | UINT64_C(0x0010000000000000)) << 10;

    int16_t const expDiff = expZ - expC;
    uint32_t sig128C[4];

    if (expDiff < 0) {
        expZ = expC;

        if ((signZ == signC) || (expDiff < -1)) {
            shiftDist -= expDiff;

            if (shiftDist) {
                sigZ = softfloat_shiftRightJam64(sigZ, static_cast<uint32_t>(shiftDist));
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

        if (0 == expDiff) {
            sigZ = static_cast<uint64_t>(sig128Z[indexWord(4, 3)]) << 32 | sig128Z[indexWord(4, 2)];
        } else {
            sig128C[indexWord(4, 3)] = sigC >> 32;
            sig128C[indexWord(4, 2)] = static_cast<uint32_t>(sigC);
            sig128C[indexWord(4, 1)] = 0;
            sig128C[indexWord(4, 0)] = 0;
            softfloat_shiftRightJam128M(sig128C, static_cast<uint8_t>(expDiff), sig128C);
        }
    }

    if (signZ == signC) {
        if (expDiff <= 0) {
            sigZ += sigC;
        } else {
            softfloat_add128M(sig128Z, sig128C, sig128Z);
            sigZ =
                static_cast<uint64_t>(sig128Z[indexWord(4, 3)]) << 32 |
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
            }

            sig128C[indexWord(4, 3)] = sigC >> 32;
            sig128C[indexWord(4, 2)] = static_cast<uint32_t>(sigC);
            sig128C[indexWord(4, 1)] = 0;
            sig128C[indexWord(4, 0)] = 0;
            softfloat_sub128M(sig128C, sig128Z, sig128Z);
        } else if (!expDiff) {
            sigZ -= sigC;

            if (0 == sigZ && !sig128Z[indexWord(4, 1)] && !sig128Z[indexWord(4, 0)]) {
                return u_as_f_64(packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
            }

            sig128Z[indexWord(4, 3)] = sigZ >> 32;
            sig128Z[indexWord(4, 2)] = static_cast<uint32_t>(sigZ);

            if (sigZ & INT64_MIN) {
                signZ = !signZ;
                softfloat_negX128M(sig128Z);
            }
        } else {
            softfloat_sub128M(sig128Z, sig128C, sig128Z);

            if (1 < expDiff) {
                sigZ = static_cast<uint64_t>(sig128Z[indexWord(4, 3)]) << 32 | sig128Z[indexWord(4, 2)];

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
        sigZ = static_cast<uint64_t>(sig128Z[indexWord(4, 3)]) << 32 | sig128Z[indexWord(4, 2)];

        if (0 == sigZ) {
            shiftDist = 64;
            sigZ = static_cast<uint64_t>(sig128Z[indexWord(4, 1)]) << 32 | sig128Z[indexWord(4, 0)];
        }

        shiftDist += count_leading_zeros(sigZ) - 1;

        if (0 != shiftDist) {
            expZ -= shiftDist;
            softfloat_shiftLeft128M(sig128Z, static_cast<uint32_t>(shiftDist), sig128Z);
            sigZ = static_cast<uint64_t>(sig128Z[indexWord(4, 3)]) << 32 | sig128Z[indexWord(4, 2)];
        }
    }

    return
        softfloat_roundPackToF64(signZ,
                                 expZ - 1,
                                 sigZ | (sig128Z[indexWord(4, 1)] || sig128Z[indexWord(4, 0)]));
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
    return mulAdd(softfloat_mulAdd_madd, f_as_u_64(a), f_as_u_64(b), f_as_u_64(c));
}
