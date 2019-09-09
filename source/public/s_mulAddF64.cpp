
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

#if (SOFTFLOAT_FAST_INT64)

float64_t
f64_mulAdd(float64_t a,
           float64_t b,
           float64_t c)
{
    using namespace softfloat::internals::fast_int64;

    if (is_sNaN(c) || is_NaN(a) || is_NaN(b)) {
        return propagate_NaN(propagate_NaN(a, b), c);
    }

    bool const signA = is_sign(a);
    bool const signB = is_sign(b);
    bool const signC = is_sign(c);
    bool signProd = signA != signB;

    if (is_inf(a) || is_inf(b)) {
        /* a or b is inf, product is inf or undefined, check other operand for zero */
        bool const is_product_undefined = is_inf(a) ? is_zero(b) : is_zero(a);

        if (is_product_undefined) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return propagate_NaN(u_as_f(defaultNaNF64UI), c);
        }

        if (is_NaN(c)) {
            return propagate_NaN(u_as_f(defaultNaNF64UI), c);
        }

        if (is_finite(c) || signProd == signC) {
            /* if summand c is finite or same sign return product as result*/
            return make_signed_inf<float64_t>(signProd);
        }

        /* summands are different sign inf, undefined sum */
        softfloat_raiseFlags(softfloat_flag_invalid);
        return propagate_NaN(u_as_f(defaultNaNF64UI), c);
    }

    if (is_NaN(c)) {
        return propagate_NaN(u_as_f(defaultNaNF64UI), c);
    }

    if (is_inf(c)) {
        return c;
    }

    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    if (is_zero(a) || is_zero(b)) {
        return is_zero(c) && signProd != signC ?
               make_signed_zero<float64_t>(softfloat_round_min == softfloat_roundingMode) :
               c;
    }

    int16_t expA = get_exp(a);
    uint64_t sigA = get_frac(a);

    if (0 == expA) {
        // subnormal or zero A

        exp16_sig64 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }


    int16_t expB = get_exp(b);
    uint64_t sigB = get_frac(b);

    if (0 == expB) {
        exp16_sig64 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int16_t expZ = expA + expB - 0x3FE;
    sigA = (sigA | UINT64_C(0x0010000000000000)) << 10;
    sigB = (sigB | UINT64_C(0x0010000000000000)) << 10;
    uint128 sig128Z = mul64To128(sigA, sigB);

    if (sig128Z.v64 < UINT64_C(0x2000000000000000)) {
        --expZ;
        sig128Z = add(sig128Z, sig128Z);
    }

    if (is_zero(c)) {
        return roundPackToF64(signProd,
                                        expZ - 1,
                                        sig128Z.v64 << 1 | (sig128Z.v0 != 0));
    }

    int16_t expC = get_exp(c);
    uint64_t sigC = get_frac(c);

    if (0 == expC) {
        exp16_sig64 const normExpSig(sigC);
        expC = normExpSig.exp;
        sigC = normExpSig.sig;
    }

    sigC = (sigC | UINT64_C(0x0010000000000000)) << 9;

    int16_t const expDiff = expZ - expC;

    uint128 sig128C;

    if (expDiff < 0) {
        expZ = expC;

        if (signProd == signC || expDiff < -1) {
            sig128Z.v64 = shiftRightJam64(sig128Z.v64, static_cast<uint32_t>(-expDiff));
        } else {
            sig128Z = shortShiftRightJam128(sig128Z, 1);
        }
    } else if (0 < expDiff) {
        sig128C = shiftRightJam128(sigC, 0u, static_cast<uint32_t>(expDiff));
    }

    /**
    @todo check case of 0 == expDiff
    */

    if (signProd == signC) {
        uint64_t sigZ;

        if (expDiff <= 0) {
            sigZ = (sigC + sig128Z.v64) | (sig128Z.v0 != 0);
        } else {
            /* 0 < expDiff */
            assert(0 < expDiff);
            sig128Z = add(sig128Z, sig128C);
            sigZ = sig128Z.v64 | (sig128Z.v0 != 0);
        }

        if (sigZ < UINT64_C(0x4000000000000000)) {
            --expZ;
            sigZ <<= 1;
        }

        return roundPackToF64(signProd, expZ, sigZ);
    }

    if (expDiff < 0) {
        signProd = signC;
        sig128Z = sub(uint128{sigC, 0}, sig128Z);
    } else if (0 == expDiff) {
        sig128Z.v64 = sig128Z.v64 - sigC;

        if (0 == (sig128Z.v64 | sig128Z.v0)) {
            return u_as_f(packToF64UI(softfloat_round_min == softfloat_roundingMode, 0, 0));
        }

        if (0 != (sig128Z.v64 & UINT64_C(0x8000000000000000))) {
            signProd = !signProd;
            sig128Z = sub(uint128{0, 0}, sig128Z);
        }
    } else {
        assert(0 < expDiff);
        sig128Z = sub(sig128Z, sig128C);
    }

    if (0 == sig128Z.v64) {
        expZ -= 64;
        sig128Z.v64 = sig128Z.v0;
        sig128Z.v0 = 0;
    }

    int16_t const shiftDist = count_leading_zeros(sig128Z.v64) - 1;
    int16_t const expZ_1 = expZ - shiftDist;

    if (shiftDist < 0) {
        uint64_t const sigZ = shortShiftRightJam64(sig128Z.v64, static_cast<uint8_t>(-shiftDist));
        return
            roundPackToF64(signProd,
                                     expZ_1,
                                     sigZ | !!(0 != sig128Z.v0));
    } else {
        uint128 const sig128Z_1 = shortShiftLeft128(sig128Z, static_cast<uint8_t>(shiftDist));
        return roundPackToF64(signProd,
                                        expZ_1,
                                        sig128Z_1.v64 | !!(0 != sig128Z_1.v0));
    }
}

#else

float64_t
f64_mulAdd(float64_t const a,
           float64_t const b,
           float64_t const c)
{
    using namespace softfloat::internals::slow_int64;

    if (is_sNaN(c) || is_NaN(a) || is_NaN(b)) {
        return propagate_NaN(propagate_NaN(a, b), c);
    }

    bool const signA = is_sign(a);
    bool const signB = is_sign(b);
    bool const signC = is_sign(c);
    bool signProd = signA != signB;

    if (is_inf(a) || is_inf(b)) {
        /* a or b is inf, product is inf or undefined, check other operand for zero */
        bool const is_product_undefined = is_inf(a) ? is_zero(b) : is_zero(a);

        if (is_product_undefined) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return propagate_NaN(u_as_f(defaultNaNF64UI), c);
        }

        if (is_NaN(c)) {
            return propagate_NaN(u_as_f(defaultNaNF64UI), c);
        }

        if (is_finite(c) || signProd == signC) {
            /* if summand c is finite or same sign return product as result*/
            return make_signed_inf<float64_t>(signProd);
        }

        /* summands are different sign inf, undefined sum */
        softfloat_raiseFlags(softfloat_flag_invalid);
        return propagate_NaN(u_as_f(defaultNaNF64UI), c);
    }

    if (is_NaN(c)) {
        return propagate_NaN(u_as_f(defaultNaNF64UI), c);
    }

    if (is_inf(c)) {
        return c;
    }

    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
    int16_t expA = get_exp(a);
    uint64_t sigA = get_frac(a);

    if (0 == expA) {
        if (0 == sigA) {
            return is_zero(c) && signProd != signC ?
                   make_signed_zero<float64_t>(softfloat_round_min == softfloat_roundingMode) :
                   c;
        }

        exp16_sig64 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    int16_t expB = get_exp(b);
    uint64_t sigB = get_frac(b);

    if (0 == expB) {
        if (0 == sigB) {
            return is_zero(c) && signProd != signC ?
                   make_signed_zero<float64_t>(softfloat_round_min == softfloat_roundingMode) :
                   c;
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

    int16_t expC = get_exp(c);
    uint64_t sigC = get_frac(c);

    if (0 == expC) {
        if (0 == sigC) {
            if (0 != shiftDist) {
                sigZ <<= 1;
            }

            return
                roundPackToF64(signProd,
                                         expZ - 1,
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

        if (signProd == signC || expDiff < -1) {
            shiftDist -= expDiff;

            if (shiftDist) {
                sigZ = shiftRightJam64(sigZ, static_cast<uint32_t>(shiftDist));
            }
        } else if (0 == shiftDist) {
            shortShiftRight128M(sig128Z, 1, sig128Z);
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

    if (signProd == signC) {
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
            sigZ = shortShiftRightJam64(sigZ, 1);
        }
    } else {
        if (expDiff < 0) {
            signProd = signC;

            if (expDiff < -1) {
                sigZ = sigC - sigZ;

                if (sig128Z[indexWord(4, 1)] || sig128Z[indexWord(4, 0)]) {
                    sigZ = (sigZ - 1) | 1;
                }

                if (!(sigZ & UINT64_C(0x4000000000000000))) {
                    --expZ;
                    sigZ <<= 1;
                }

                return roundPackToF64(signProd, expZ - 1, sigZ);
            }

            sig128C[indexWord(4, 3)] = sigC >> 32;
            sig128C[indexWord(4, 2)] = static_cast<uint32_t>(sigC);
            sig128C[indexWord(4, 1)] = 0;
            sig128C[indexWord(4, 0)] = 0;
            softfloat_sub128M(sig128C, sig128Z, sig128Z);
        } else if (!expDiff) {
            sigZ -= sigC;

            if (0 == sigZ && !sig128Z[indexWord(4, 1)] && !sig128Z[indexWord(4, 0)]) {
                return u_as_f(packToF64UI(softfloat_roundingMode == softfloat_round_min, 0, 0));
            }

            sig128Z[indexWord(4, 3)] = sigZ >> 32;
            sig128Z[indexWord(4, 2)] = static_cast<uint32_t>(sigZ);

            if (sigZ & INT64_MIN) {
                signProd = !signProd;
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
                    roundPackToF64(signProd, expZ - 1,
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
        roundPackToF64(signProd,
                                 expZ - 1,
                                 sigZ | (sig128Z[indexWord(4, 1)] || sig128Z[indexWord(4, 0)]));
}

#endif
