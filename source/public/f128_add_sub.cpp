
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
All rights reserved.

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

#include <utility>

namespace {
static float128_t
add_magnitudes(uint64_t const uiA64,
               uint64_t const uiA0,
               uint64_t const uiB64,
               uint64_t const uiB0,
               bool const signZ)
{
    using namespace softfloat::internals::fast_int64;

    int32_t const expA = exp_F128_UI64(uiA64);
    uint128 sigA{frac_F128_UI64(uiA64), uiA0};
    int32_t const expB = exp_F128_UI64(uiB64);
    uint128 sigB{frac_F128_UI64(uiB64), uiB0};
    int32_t expDiff = expA - expB;

    if (0 == expDiff) {
        if (0x7FFF == expA) {
            if (0 != (sigA.v64 | sigA.v0 | sigB.v64 | sigB.v0)) {
                return float128_t(propagate_NaN(uiA64, uiA0, uiB64, uiB0));
            }

            return static_cast<float128_t>(uint128{uiA64, uiA0});
        }

        uint128 const sigZ = add(sigA, sigB);

        if (0 == expA) {
            return static_cast<float128_t>(uint128{pack_to_F128_UI64(signZ, 0, sigZ.v64), sigZ.v0});
        }

        uint128_extra const sig128Extra = short_shift_right_jam_128Extra(uint128{sigZ.v64 | UINT64_C(0x0002000000000000), sigZ.v0}, 0, 1);
        return round_pack_to_F128(signZ, expA, sig128Extra.v.v64, sig128Extra.v.v0, sig128Extra.extra);
    }

    if (expDiff < 0) {
        if (0x7FFF == expB) {
            if (sigB.v64 | sigB.v0) {
                return float128_t(propagate_NaN(uiA64, uiA0, uiB64, uiB0));
            }

            return static_cast<float128_t>(uint128{pack_to_F128_UI64(signZ, 0x7FFF, 0), 0});
        }

        if (0 != expA) {
            sigA.v64 |= UINT64_C(0x0001000000000000);
        } else if (0 == expDiff + 1) {
            uint128 const sigZ = add(uint128{sigA.v64 | UINT64_C(0x0001000000000000), sigA.v0}, sigB);

            if (sigZ.v64 < UINT64_C(0x0002000000000000)) {
                return round_pack_to_F128(signZ, expB - 1, sigZ.v64, sigZ.v0, 0);
            }

            uint128_extra const sig128Extra = short_shift_right_jam_128Extra(sigZ, 0, 1);
            return round_pack_to_F128(signZ, expB, sig128Extra.v.v64, sig128Extra.v.v0, sig128Extra.extra);
        }

        uint128_extra const sig128Extra = shift_right_jam_128Extra(sigA, 0, static_cast<uint32_t>(-(expDiff + 1)));
        uint128 const sigZ = add(uint128{sig128Extra.v.v64 | UINT64_C(0x0001000000000000), sig128Extra.v.v0}, sigB);

        if (sigZ.v64 < UINT64_C(0x0002000000000000)) {
            return round_pack_to_F128(signZ, expB - 1, sigZ.v64, sigZ.v0, sig128Extra.extra);
        }

        uint128_extra const sig128Extra_1 = short_shift_right_jam_128Extra(sigZ, sig128Extra.extra, 1);
        return round_pack_to_F128(signZ, expB, sig128Extra_1.v.v64, sig128Extra_1.v.v0, sig128Extra_1.extra);
    }

    if (0x7FFF == expA) {
        if (0 != (sigA.v64 | sigA.v0)) {
            return float128_t(propagate_NaN(uiA64, uiA0, uiB64, uiB0));
        }

        return static_cast<float128_t>(uint128{uiA64, uiA0});
    }

    if (0 != expB) {
        sigB.v64 |= UINT64_C(0x0001000000000000);
    } else if (1 == expDiff) {
        uint128 const sigZ = add(uint128{sigA.v64 | UINT64_C(0x0001000000000000), sigA.v0}, sigB);

        if (sigZ.v64 < UINT64_C(0x0002000000000000)) {
            return round_pack_to_F128(signZ, expA - 1, sigZ.v64, sigZ.v0, 0);
        }

        uint128_extra const sig128Extra = short_shift_right_jam_128Extra(sigZ, 0, 1);
        return round_pack_to_F128(signZ, expA, sig128Extra.v.v64, sig128Extra.v.v0, sig128Extra.extra);
    }

    uint128_extra const sig128Extra = shift_right_jam_128Extra(sigB, 0, static_cast<uint32_t>(expDiff - 1));
    uint128 const sigZ = add(uint128{sigA.v64 | UINT64_C(0x0001000000000000), sigA.v0}, sig128Extra.v);

    if (sigZ.v64 < UINT64_C(0x0002000000000000)) {
        return round_pack_to_F128(signZ, expA - 1, sigZ.v64, sigZ.v0, sig128Extra.extra);
    }

    uint128_extra const sig128Extra_1 = short_shift_right_jam_128Extra(sigZ, sig128Extra.extra, 1);
    return round_pack_to_F128(signZ, expA, sig128Extra_1.v.v64, sig128Extra_1.v.v0, sig128Extra_1.extra);
}

static float128_t
sub_magnitudes(uint64_t const uiA64,
               uint64_t const uiA0,
               uint64_t const uiB64,
               uint64_t const uiB0,
               bool const signZ)
{
    using namespace softfloat::internals::fast_int64;

    int32_t const expA = exp_F128_UI64(uiA64);
    int32_t const expB = exp_F128_UI64(uiB64);
    uint128 const sigA = short_shift_left_128(uint128{frac_F128_UI64(uiA64), uiA0}, 4);
    uint128 const sigB = short_shift_left_128(uint128{frac_F128_UI64(uiB64), uiB0}, 4);
    int32_t const expDiff = expA - expB;

    if (!(0 < expDiff)) {
        if (!(expDiff < 0)) {

            if (0x7FFF == expA) {
                if (0 != (sigA.v64 | sigA.v0 | sigB.v64 | sigB.v0)) {
                    return float128_t(propagate_NaN(uiA64, uiA0, uiB64, uiB0));
                }

                softfloat_raiseFlags(softfloat_flag_invalid);
                return float128_t(uint128{defaultNaNF128UI64, defaultNaNF128UI0});
            }

            int32_t const expZ = 0 == expA ? 1 : expA;

            if (sigB.v64 < sigA.v64) {
                uint128 const sigZ = sub(sigA, sigB);
                return norm_round_pack_to_F128(signZ, expZ - 5, sigZ.v64, sigZ.v0);
            }

            if (sigA.v64 < sigB.v64) {
                uint128 const sigZ = sub(sigB, sigA);
                return norm_round_pack_to_F128(!signZ, expZ - 5, sigZ.v64, sigZ.v0);
            }

            if (sigB.v0 < sigA.v0) {
                uint128 const sigZ = sub(sigA, sigB);
                return norm_round_pack_to_F128(signZ, expZ - 5, sigZ.v64, sigZ.v0);
            }

            if (sigA.v0 < sigB.v0) {
                uint128 const sigZ = sub(sigB, sigA);
                return norm_round_pack_to_F128(!signZ, expZ - 5, sigZ.v64, sigZ.v0);
            }

            softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
            return float128_t(uint128{pack_to_F128_UI64(softfloat_roundingMode == softfloat_round_min, 0, 0), 0});
        }

        if (0x7FFF == expB) {
            if (0 != (sigB.v64 | sigB.v0)) {
                return float128_t(propagate_NaN(uiA64, uiA0, uiB64, uiB0));
            }

            return float128_t(uint128{pack_to_F128_UI64(!signZ, 0x7FFF, 0), 0});
        }

        if (0 != expA) {
            auto const sigA_1 = shift_right_jam_128(sigA.v64 | UINT64_C(0x0010000000000000), sigA.v0, static_cast<uint32_t>(-expDiff));
            uint128 const sigZ = sub(uint128{sigB.v64 | UINT64_C(0x0010000000000000), sigB.v0}, sigA_1);
            return norm_round_pack_to_F128(!signZ, expB - 5, sigZ.v64, sigZ.v0);
        }

        auto const expDiff_1 = expDiff + 1;

        if (0 == expDiff_1) {
            uint128 const sigZ = sub(uint128{sigB.v64 | UINT64_C(0x0010000000000000), sigB.v0}, sigA);
            return norm_round_pack_to_F128(!signZ, expB - 5, sigZ.v64, sigZ.v0);
        }

        auto const sigA_1 = shift_right_jam_128(sigA, static_cast<uint32_t>(-expDiff_1));
        uint128 const sigZ = sub(uint128{sigB.v64 | UINT64_C(0x0010000000000000), sigB.v0}, sigA_1);
        return norm_round_pack_to_F128(!signZ, expB - 5, sigZ.v64, sigZ.v0);
    }

    if (0x7FFF == expA) {
        if (0 != (sigA.v64 | sigA.v0)) {
            return float128_t(propagate_NaN(uiA64, uiA0, uiB64, uiB0));
        }

        return float128_t(uint128{uiA64, uiA0});
    }

    if (0 != expB) {
        auto const sigB_1 = shift_right_jam_128(sigB.v64 | UINT64_C(0x0010000000000000), sigB.v0, static_cast<uint32_t>(expDiff));
        uint128 const sigZ = sub(uint128{sigA.v64 | UINT64_C(0x0010000000000000), sigA.v0}, sigB_1);
        return norm_round_pack_to_F128(signZ, expA - 5, sigZ.v64, sigZ.v0);
    }

    auto const expDiff_1 = expDiff - 1;

    if (0 == expDiff_1) {
        uint128 const sigZ = sub(uint128{sigA.v64 | UINT64_C(0x0010000000000000), sigA.v0}, sigB);
        return norm_round_pack_to_F128(signZ, expA - 5, sigZ.v64, sigZ.v0);
    }

    auto const sigB_1 = shift_right_jam_128(sigB, static_cast<uint32_t>(expDiff_1));
    uint128 const sigZ = sub(uint128{sigA.v64 | UINT64_C(0x0010000000000000), sigA.v0}, sigB_1);
    return norm_round_pack_to_F128(signZ, expA - 5, sigZ.v64, sigZ.v0);
}
#if !(SOFTFLOAT_FAST_INT64)
static void
add_M_F128(uint32_t const* aWPtr,
           uint32_t const* bWPtr,
           uint32_t* const zWPtr,
           bool negateB)
{
    using namespace softfloat::internals::slow_int64;

    uint32_t uiA96 = aWPtr[index_word_hi(4)];
    int32_t expA = exp_F128_UI96(uiA96);
    uint32_t uiB96 = bWPtr[index_word_hi(4)];
    int32_t expB = exp_F128_UI96(uiB96);

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (try_propagate_NaN_M_F128(aWPtr, bWPtr, zWPtr)) {
            return;
        }

        uint32_t uiZ96 = uiA96;

        if (expB == 0x7FFF) {
            uiZ96 = uiB96 ^ pack_to_F128_UI96(negateB, 0, 0);

            if (0x7FFF == expA && uiZ96 != uiA96) {
                invalid_M_F128(zWPtr);
                return;
            }
        }

        zWPtr[index_word_hi(4)] = uiZ96;
        zWPtr[index_word(4, 2)] = 0;
        zWPtr[index_word(4, 1)] = 0;
        zWPtr[index_word(4, 0)] = 0;
        return;
    }

    bool signZ = is_sign(uiA96);
    bool const signB = is_sign(uiB96) != negateB;
    negateB = signZ != signB;

    if (static_cast<uint32_t>(uiA96 << 1) < static_cast<uint32_t>(uiB96 << 1)) {
        signZ = signB;
        expA = expB;
        expB = exp_F128_UI96(uiA96);
        std::swap(aWPtr, bWPtr);
        uiA96 = uiB96;
        uiB96 = bWPtr[index_word_hi(4)];
    }

    uint32_t sig96A = frac_F128_UI96(uiA96);
    uint32_t sig96B = frac_F128_UI96(uiB96);

    if (expA) {
        --expA;
        sig96A |= 0x00010000;

        if (expB) {
            --expB;
            sig96B |= 0x00010000;
        }
    }

    auto const addCarryMRoutinePtr = negateB ? add_compl_carry_M : add_carry_M;
    int32_t const expDiff = expA - expB;

    bool carry;
    uint32_t wordSigZ;
    uint32_t extSigZ[5];

    if (expDiff) {
        extSigZ[index_word_hi(5)] = sig96B;
        extSigZ[index_word(5, 3)] = bWPtr[index_word(4, 2)];
        extSigZ[index_word(5, 2)] = bWPtr[index_word(4, 1)];
        extSigZ[index_word(5, 1)] = bWPtr[index_word(4, 0)];
        extSigZ[index_word(5, 0)] = 0;
        shift_right_jam_M_160(extSigZ, static_cast<uint8_t>(expDiff), extSigZ);
        sig96B = extSigZ[index_word_hi(5)];
        carry = 0;

        if (negateB) {
            sig96B = ~sig96B;
            wordSigZ = extSigZ[index_word_lo(5)];
            extSigZ[index_word_lo(5)] = static_cast<uint32_t>(-static_cast<int32_t>(wordSigZ));
            carry = !wordSigZ;
        }

        carry = (*addCarryMRoutinePtr)(3,
                                       &aWPtr[index_multiword_lo(4, 3)],
                                       &extSigZ[index_multiword(5, 3, 1)],
                                       carry,
                                       &extSigZ[index_multiword(5, 3, 1)]);
        wordSigZ = sig96A + sig96B + !!(carry);
    } else {
        extSigZ[index_word_lo(5)] = 0;
        carry = (*addCarryMRoutinePtr)(3,
                                       &aWPtr[index_multiword_lo(4, 3)],
                                       &bWPtr[index_multiword_lo(4, 3)],
                                       negateB,
                                       &extSigZ[index_multiword(5, 3, 1)]);

        if (negateB) {
            wordSigZ = sig96A + ~sig96B + !!carry;

            if (wordSigZ & 0x80000000) {
                signZ = !signZ;
                carry = add_compl_carry_M_96(&bWPtr[index_multiword_lo(4, 3)],
                                             &aWPtr[index_multiword_lo(4, 3)],
                                             1,
                                             &extSigZ[index_multiword(5, 3, 1)]);
                wordSigZ = sig96B + ~sig96A + !!carry;
            } else {
                if (0 == wordSigZ && 0 == extSigZ[index_word(5, 3)] && 0 == (extSigZ[index_word(5, 2)] | extSigZ[index_word(5, 1)] | extSigZ[index_word(5, 0)])) {
                    softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
                    signZ = softfloat_round_min == softfloat_roundingMode;
                    zWPtr[index_word_hi(4)] = pack_to_F128_UI96(signZ, 0, 0);
                    zWPtr[index_word(4, 2)] = 0;
                    zWPtr[index_word(4, 1)] = 0;
                    zWPtr[index_word(4, 0)] = 0;
                    return;
                }
            }
        } else {
            wordSigZ = sig96A + sig96B + !!carry;
        }
    }

    extSigZ[index_word_hi(5)] = wordSigZ;

    if (0x00010000 <= wordSigZ) {
        if (0x00020000 <= wordSigZ) {
            ++expA;
            short_shift_right_jam_M_160(extSigZ, 1, extSigZ);
        }

        round_pack_to_M_F128(signZ, expA, extSigZ, zWPtr);
    } else {
        norm_round_pack_to_M_F128(signZ, expA, extSigZ, zWPtr);
    }

}
#endif
}  // namespace

float128_t
f128_add(float128_t const a,
         float128_t const b)
{
    using namespace softfloat::internals::fast_int64;

    uint128 const uA{a};
    uint128 const uB{b};
    bool const signA = is_sign(uA.v64);
    bool const signB = is_sign(uB.v64);
    return
        signA == signB ?
        add_magnitudes(uA.v64, uA.v0, uB.v64, uB.v0, signA) :
        sub_magnitudes(uA.v64, uA.v0, uB.v64, uB.v0, signA);
}

float128_t
f128_sub(float128_t const a,
         float128_t const b)
{
    using namespace softfloat::internals::fast_int64;

    uint128 const uA{a};
    bool const signA = is_sign(uA.v64);
    uint128 const uB{b};
    bool const signB = is_sign(uB.v64);
    return
        signA == signB ?
        sub_magnitudes(uA.v64, uA.v0, uB.v64, uB.v0, signA) :
        add_magnitudes(uA.v64, uA.v0, uB.v64, uB.v0, signA);
}

void
f128M_add(float128_t const *const aPtr,
          float128_t const *const bPtr,
          float128_t* const zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    using namespace softfloat::internals::fast_int64;

    uint64_t const* const aWPtr = reinterpret_cast<uint64_t const*>(aPtr);
    uint64_t const* const bWPtr = reinterpret_cast<uint64_t const*>(bPtr);
    uint64_t const uiA64 = aWPtr[index_word(2, 1)];
    uint64_t const uiA0 = aWPtr[index_word(2, 0)];
    bool const signA = is_sign(uiA64);

    uint64_t const uiB64 = bWPtr[index_word(2, 1)];
    uint64_t const uiB0 = bWPtr[index_word(2, 0)];
    bool const signB = is_sign(uiB64);

    *zPtr =
        signA == signB ?
        add_magnitudes(uiA64, uiA0, uiB64, uiB0, signA) :
        sub_magnitudes(uiA64, uiA0, uiB64, uiB0, signA);
#else
    using namespace softfloat::internals::slow_int64;
    add_M_F128(reinterpret_cast<const uint32_t*>(aPtr),
               reinterpret_cast<const uint32_t*>(bPtr),
               reinterpret_cast<uint32_t*>(zPtr),
               false);
#endif
}

void
f128M_sub(float128_t const* const aPtr,
          float128_t const* const bPtr,
          float128_t* const  zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    using namespace softfloat::internals::fast_int64;

    auto const aWPtr = reinterpret_cast<uint64_t const*>(aPtr);
    auto const bWPtr = reinterpret_cast<uint64_t const*>(bPtr);
    uint64_t const uiA64 = aWPtr[index_word(2, 1)];
    uint64_t const uiA0 = aWPtr[index_word(2, 0)];
    uint64_t const uiB64 = bWPtr[index_word(2, 1)];
    uint64_t const uiB0 = bWPtr[index_word(2, 0)];
    bool const signA = is_sign(uiA64);
    bool const signB = is_sign(uiB64);

    *zPtr =
        signA == signB ?
        sub_magnitudes(uiA64, uiA0, uiB64, uiB0, signA) :
        add_magnitudes(uiA64, uiA0, uiB64, uiB0, signA);

#else
    using namespace softfloat::internals::slow_int64;

    add_M_F128(reinterpret_cast<uint32_t const *>(aPtr),
               reinterpret_cast<uint32_t const *>(bPtr),
               reinterpret_cast<uint32_t*>(zPtr),
               true);
#endif
}
