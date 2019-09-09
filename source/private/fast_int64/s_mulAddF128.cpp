
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

#include <cassert>
#include <cstring>

namespace softfloat {
namespace internals {
namespace fast_int64 {
namespace {

static
void
short_shift_right_jam_M(size_t const size_words,
                        uint64_t const* const aPtr,
                        uint8_t const dist,
                        uint64_t* const zPtr)
{
    uint8_t const uNegDist = 63u & -static_cast<int8_t>(dist);
    auto index = index_word_lo(size_words);
    auto const lastIndex = index_word_hi(size_words);
    uint64_t wordA = aPtr[index];
    uint64_t partWordZ = wordA >> dist;

    if (partWordZ << dist != wordA) {
        partWordZ |= 1;
    }

    while (index != lastIndex) {
        wordA = aPtr[index + wordIncr];
        zPtr[index] = wordA << uNegDist | partWordZ;
        index += wordIncr;
        partWordZ = wordA >> dist;
    }

    zPtr[index] = partWordZ;

}

/**
Shifts the 256-bit unsigned integer pointed to by `aPtr' right by the number
of bits given in `dist', which must not be zero.  If any nonzero bits are
shifted off, they are "jammed" into the least-significant bit of the shifted
value by setting the least-significant bit to 1.  This shifted-and-jammed
value is stored at the location pointed to by `zPtr'.  Each of `aPtr' and
`zPtr' points to an array of four 64-bit elements that concatenate in the
platform's normal endian order to form a 256-bit integer.
The value of `dist' can be arbitrarily large.  In particular, if `dist'
is greater than 256, the stored result will be either 0 or 1, depending on
whether the original 256 bits are all zeros.
*/
static void
shift_right_jam_M_256(uint64_t const* aPtr,
                            uint32_t const dist,
                            uint64_t* const zPtr)
{
    assert(aPtr);
    assert(zPtr);

    uint64_t wordJam = 0;
    uint64_t* ptr = nullptr;
    uint32_t wordDist = dist / (CHAR_BIT * sizeof(*aPtr));

    if (0 != wordDist) {
        if (4 < wordDist) {
            wordDist = 4;
        }

        assert(wordDist <= 4);

        uint64_t const* ptr1 = aPtr + index_multiword_lo(4, wordDist);

        for (uint32_t i = wordDist; 0 != i; --i) {
            wordJam = *ptr1++;

            if (0 != wordJam) {
                break;
            }
        }

        ptr = zPtr;
    }

    assert(wordDist <= 4);

    if (wordDist < 4) {
        aPtr += index_multiword_hi_but(4, wordDist);
        uint8_t const innerDist = dist & 63;

        if (0 != innerDist) {
            short_shift_right_jam_M(4u - wordDist,
                                    aPtr,
                                    innerDist,
                                    zPtr + index_multiword_lo_but(4, wordDist));

            if (0 == wordDist) {
                if (wordJam) {
                    zPtr[index_word_lo(4)] |= 1;
                }

                return;
            }
        } else {
            aPtr += index_word_lo(4 - wordDist);
            ptr = zPtr + index_word_lo(4);

            /**
            @todo memmove
            */
            for (auto i = 4 - wordDist; i; --i, aPtr += wordIncr, ptr += wordIncr) {
                *ptr = *aPtr;
            }
        }

        ptr = zPtr + index_multiword_hi(4, wordDist);
    }

    assert(ptr);
    memset(ptr, 0, wordDist * sizeof(*ptr));

    if (wordJam) {
        zPtr[index_word_lo(4)] |= 1;
    }
}

}  // namespace

/**
@todo Merge fast int64 and non fast int 64
*/
float128_t
mul_add_F128(Mul_add_operations const op,
             uint64_t const& uiA64,
             uint64_t const& uiA0,
             uint64_t const& uiB64,
             uint64_t const& uiB0,
             uint64_t const& uiC64,
             uint64_t const& uiC0)
{
    bool const signA = is_sign(uiA64);
    int32_t expA = exp_F128_UI64(uiA64);
    uint128 sigA{frac_F128_UI64(uiA64), uiA0};
    bool const signB = is_sign(uiB64);
    int32_t expB = exp_F128_UI64(uiB64);
    uint128 sigB{frac_F128_UI64(uiB64), uiB0};
    bool const signC = is_sign(uiC64) != (softfloat_mulAdd_subC == op);
    int32_t expC = exp_F128_UI64(uiC64);
    uint128 sigC{frac_F128_UI64(uiC64), uiC0};
    bool signZ = (signA != signB) != (softfloat_mulAdd_subProd == op);

    if (0x7FFF == expA) {
        if (0 != (sigA.v64 | sigA.v0) || (expB == 0x7FFF && 0 != (sigB.v64 | sigB.v0))) {
            auto const uiZ_1 = propagate_NaN(uiA64, uiA0, uiB64, uiB0);
            return float128_t(propagate_NaN(uiZ_1.v64, uiZ_1.v0, uiC64, uiC0));
        }

        if (0 != (expB | sigB.v64 | sigB.v0)) {
            uint128 const uiZ_1 = uint128{pack_to_F128_UI64(signZ, 0x7FFF, 0), 0};

            if (0x7FFF != expC) {
                return static_cast<float128_t>(uiZ_1);
            }

            if (0 != (sigC.v64 | sigC.v0)) {
                return float128_t(propagate_NaN(uiZ_1.v64, uiZ_1.v0, uiC64, uiC0));
            }

            if (signZ == signC) {
                return static_cast<float128_t>(uiZ_1);
            }
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return static_cast<float128_t>(propagate_NaN(defaultNaNF128UI64, defaultNaNF128UI0, uiC64, uiC0));
    }

    if (0x7FFF == expB) {
        if (0 != (sigB.v64 | sigB.v0)) {
            auto const uiZ_1 = propagate_NaN(uiA64, uiA0, uiB64, uiB0);
            return float128_t(propagate_NaN(uiZ_1.v64, uiZ_1.v0, uiC64, uiC0));
        }

        if (0 != (expA | sigA.v64 | sigA.v0)) {
            auto const uiZ_1 = uint128{pack_to_F128_UI64(signZ, 0x7FFF, 0), 0};

            if (0x7FFF != expC) {
                return static_cast<float128_t>(uiZ_1);
            }

            if (0 != (sigC.v64 | sigC.v0)) {
                return float128_t(propagate_NaN(uiZ_1.v64, uiZ_1.v0, uiC64, uiC0));
            }

            if (signZ == signC) {
                return float128_t(uiZ_1);
            }
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return float128_t(propagate_NaN(defaultNaNF128UI64, defaultNaNF128UI0, uiC64, uiC0));
    }

    if (0x7FFF == expC) {
        if (0 != (sigC.v64 | sigC.v0)) {
            return float128_t(propagate_NaN(0, 0, uiC64, uiC0));
        }

        return static_cast<float128_t>(uint128{uiC64, uiC0});
    }

    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();

    if (0 == expA) {
        if (0 == (sigA.v64 | sigA.v0)) {
            if (0 == (expC | sigC.v64 | sigC.v0) && signZ != signC) {
                return static_cast<float128_t>(uint128{pack_to_F128_UI64((softfloat_roundingMode == softfloat_round_min), 0, 0), 0});
            }

            return static_cast<float128_t>(uint128{uiC64, uiC0});
        }

        exp32_sig128 const normExpSig = norm_subnormal_F128Sig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (0 == expB) {
        if (0 == (sigB.v64 | sigB.v0)) {
            if (0 == (expC | sigC.v64 | sigC.v0) && signZ != signC) {
                return static_cast<float128_t>(uint128{pack_to_F128_UI64((softfloat_roundingMode == softfloat_round_min), 0, 0), 0});
            }

            return static_cast<float128_t>(uint128{uiC64, uiC0});
        }

        exp32_sig128 const normExpSig = norm_subnormal_F128Sig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int32_t expZ = expA + expB - 0x3FFE;
    sigA.v64 |= UINT64_C(0x0001000000000000);
    sigB.v64 |= UINT64_C(0x0001000000000000);
    sigA = short_shift_left_128(sigA, 8);
    sigB = short_shift_left_128(sigB, 15);

    uint64_t sig256Z[4];
    mul_M_128_to_256(sigA, sigB, sig256Z);

    uint128 sigZ{sig256Z[index_word(4, 3)], sig256Z[index_word(4, 2)]};

    int32_t shiftDist_1 = 0;

    if (0 == (sigZ.v64 & UINT64_C(0x0100000000000000))) {
        --expZ;
        shiftDist_1 = -1;
    }

    if (0 == expC) {
        if (0 == (sigC.v64 | sigC.v0)) {
            shiftDist_1 += 8;
            auto const sigZExtra_2 = sig256Z[index_word(4, 1)] | sig256Z[index_word(4, 0)];
            auto const sigZExtra_1 = static_cast<uint64_t>(sigZ.v0 << (64 - shiftDist_1)) | (sigZExtra_2 != 0);
            auto const sigZ_1 = short_shift_right_128(sigZ, static_cast<uint8_t>(shiftDist_1));
            return round_pack_to_F128(signZ, expZ - 1, sigZ_1.v64, sigZ_1.v0, sigZExtra_1);
        }

        exp32_sig128 const normExpSig = norm_subnormal_F128Sig(sigC);
        expC = normExpSig.exp;
        sigC = normExpSig.sig;
    }

    sigC.v64 |= UINT64_C(0x0001000000000000);
    sigC = short_shift_left_128(sigC, 8);

    int32_t const expDiff = expZ - expC;
    uint64_t sig256C[4];

    if (expDiff < 0) {
        expZ = expC;

        if (signZ == signC || expDiff < -1) {
            shiftDist_1 -= expDiff;

            if (0 != shiftDist_1) {
                sigZ = shift_right_jam_128(sigZ, static_cast<uint32_t>(shiftDist_1));
            }
        } else {
            if (0 == shiftDist_1) {
                uint128 const x128 = short_shift_right_128(uint128{sig256Z[index_word(4, 1)], sig256Z[index_word(4, 0)]}, 1);
                sig256Z[index_word(4, 1)] = (sigZ.v0 << 63) | x128.v64;
                sig256Z[index_word(4, 0)] = x128.v0;
                sigZ = short_shift_right_128(sigZ, 1);
                sig256Z[index_word(4, 3)] = sigZ.v64;
                sig256Z[index_word(4, 2)] = sigZ.v0;
            }
        }
    } else {
        if (0 != shiftDist_1) {
            add_M_256(sig256Z, sig256Z, sig256Z);
        }

        if (0 == expDiff) {
            sigZ = uint128{sig256Z[index_word(4, 3)], sig256Z[index_word(4, 2)]};
        } else {
            sig256C[index_word(4, 3)] = sigC.v64;
            sig256C[index_word(4, 2)] = sigC.v0;
            sig256C[index_word(4, 1)] = 0;
            sig256C[index_word(4, 0)] = 0;
            shift_right_jam_M_256(sig256C, static_cast<uint32_t>(expDiff), sig256C);
        }
    }

    {
        int32_t shiftDist_2 = 8;

        if (signZ == signC) {
            if (expDiff <= 0) {
                sigZ = add(sigC, sigZ);
            } else {
                add_M_256(sig256Z, sig256C, sig256Z);
                sigZ = uint128{sig256Z[index_word(4, 3)], sig256Z[index_word(4, 2)]};
            }

            if (sigZ.v64 & UINT64_C(0x0200000000000000)) {
                ++expZ;
                shiftDist_2 = 9;
            }

            auto const sigZExtra_2 = sig256Z[index_word(4, 1)] | sig256Z[index_word(4, 0)];
            auto const sigZExtra_1 = static_cast<uint64_t>(sigZ.v0 << (64 - shiftDist_2)) | (sigZExtra_2 != 0);
            auto const sigZ_1 = short_shift_right_128(sigZ, static_cast<uint8_t>(shiftDist_2));
            return round_pack_to_F128(signZ, expZ - 1, sigZ_1.v64, sigZ_1.v0, sigZExtra_1);
        }

        {
            uint64_t sigZExtra;

            if (expDiff < 0) {
                signZ = signC;

                if (expDiff < -1) {
                    sigZ = sub(sigC, sigZ);
                    sigZExtra = sig256Z[index_word(4, 1)] | sig256Z[index_word(4, 0)];

                    if (0 != sigZExtra) {
                        sigZ = sub(sigZ, uint128{0, 1});
                    }

                    if (0 == (sigZ.v64 & UINT64_C(0x0100000000000000))) {
                        --expZ;
                        shiftDist_2 = 7;
                    }

                    auto const sigZExtra_1 = static_cast<uint64_t>(sigZ.v0 << (64 - shiftDist_2)) | !!(0 != sigZExtra);
                    auto const sigZ_1 = short_shift_right_128(sigZ, static_cast<uint8_t>(shiftDist_2));
                    return round_pack_to_F128(signZ, expZ - 1, sigZ_1.v64, sigZ_1.v0, sigZExtra_1);
                }

                sig256C[index_word(4, 3)] = sigC.v64;
                sig256C[index_word(4, 2)] = sigC.v0;
                sig256C[index_word(4, 1)] = 0;
                sig256C[index_word(4, 0)] = 0;
                sub_M_256(sig256C, sig256Z, sig256Z);
            } else if (0 == expDiff) {
                sigZ = sub(sigZ, sigC);

                if (0 == (sigZ.v64 | sigZ.v0) && 0 == sig256Z[index_word(4, 1)] && 0 == sig256Z[index_word(4, 0)]) {
                    return static_cast<float128_t>(uint128{pack_to_F128_UI64(softfloat_round_min == softfloat_roundingMode, 0, 0), 0});
                }

                sig256Z[index_word(4, 3)] = sigZ.v64;
                sig256Z[index_word(4, 2)] = sigZ.v0;

                if (0 != (sigZ.v64 & UINT64_C(0x8000000000000000))) {
                    signZ = !signZ;
                    uint64_t const zero256[4] = {0, 0, 0, 0};
                    sub_M_256(zero256, sig256Z, sig256Z);
                }
            } else {
                sub_M_256(sig256Z, sig256C, sig256Z);

                if (1 < expDiff) {
                    sigZ.v64 = sig256Z[index_word(4, 3)];
                    sigZ.v0 = sig256Z[index_word(4, 2)];

                    if (0 == (sigZ.v64 & UINT64_C(0x0100000000000000))) {
                        --expZ;
                        shiftDist_2 = 7;
                    }

                    auto const sigZExtra_2 = sig256Z[index_word(4, 1)] | sig256Z[index_word(4, 0)];
                    auto const sigZExtra_1 = static_cast<uint64_t>(sigZ.v0 << (64 - shiftDist_2)) | (sigZExtra_2 != 0);
                    auto const sigZ_1 = short_shift_right_128(sigZ, static_cast<uint8_t>(shiftDist_2));
                    return round_pack_to_F128(signZ, expZ - 1, sigZ_1.v64, sigZ_1.v0, sigZExtra_1);
                }
            }
        }

        sigZ = uint128{sig256Z[index_word(4, 3)], sig256Z[index_word(4, 2)]};

        {
            uint64_t sigZExtra = sig256Z[index_word(4, 1)];
            uint64_t const sig256Z0 = sig256Z[index_word(4, 0)];

            if (0 != sigZ.v64) {
                if (0 != sig256Z0) {
                    sigZExtra |= 1;
                }
            } else {
                expZ -= 64;
                sigZ = uint128{sigZ.v0, sigZExtra};
                sigZExtra = sig256Z0;

                if (0 == sigZ.v64) {
                    expZ -= 64;
                    sigZ = uint128{sigZ.v0, sigZExtra};
                    sigZExtra = 0;

                    if (0 == sigZ.v64) {
                        expZ -= 64;
                        sigZ = uint128{sigZ.v0, 0};
                    }
                }
            }

            int32_t const shiftDist_2_3 = count_leading_zeros(sigZ.v64);
            expZ += 7 - shiftDist_2_3;
            int32_t const shiftDist_2_2 = 15 - shiftDist_2_3;

            if (0 < shiftDist_2_2) {
                auto const sigZExtra_1 = static_cast<uint64_t>(sigZ.v0 << (64 - shiftDist_2_2)) | !!(0 != sigZExtra);
                auto const sigZ_1 = short_shift_right_128(sigZ, static_cast<uint8_t>(shiftDist_2_2));
                return round_pack_to_F128(signZ, expZ - 1, sigZ_1.v64, sigZ_1.v0, sigZExtra_1);
            }

            if (0 != shiftDist_2_2) {
                auto const shiftDist_2_1 = -shiftDist_2_2;
                sigZ = short_shift_left_128(sigZ, static_cast<uint8_t>(shiftDist_2_1));
                uint128 const x128 = short_shift_left_128(uint128{0, sigZExtra}, static_cast<uint8_t>(shiftDist_2_1));
                sigZ.v0 |= x128.v64;
                sigZExtra = x128.v0;
            }

            return round_pack_to_F128(signZ, expZ - 1, sigZ.v64, sigZ.v0, sigZExtra);
        }
    }
}

}  // namespace fast_int64
}  // namespace internals
}  // namespace softfloat
