
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
All rights reserved.
*/
/*
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

#include <algorithm>

namespace {

#if !(SOFTFLOAT_FAST_INT64)
static void
add_M_extF80(extFloat80M const* aSPtr,
             extFloat80M const* bSPtr,
             extFloat80M* const zSPtr,
             bool negateB)
{
    using namespace softfloat::internals::slow_int64;
    typedef void(*round_packer_ptr)(bool, int32_t, uint32_t*, uint8_t const&, extFloat80M*);

    int32_t expA = get_exp(*aSPtr);
    int32_t expB = get_exp(*bSPtr);

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (try_propagate_NaN_M_extF80(aSPtr, bSPtr, zSPtr)) {
            return;
        }

        uint16_t uiZ64 = aSPtr->signExp;

        if (0x7FFF == expB) {
            uiZ64 = static_cast<uint16_t>(bSPtr->signExp ^ pack_to_extF80_UI64(negateB, 0));

            if (0x7FFF == expA && uiZ64 != aSPtr->signExp) {
                invalid_M_extF80(zSPtr);
                return;
            }
        }

        zSPtr->signExp = uiZ64;
        zSPtr->signif = static_cast<uint64_t>(INT64_MIN);
        return;
    }

    bool signZ = is_sign(*aSPtr);
    bool const signB = is_sign(*bSPtr) != negateB;
    negateB = signZ != signB;

    if (expA < expB) {
        signZ = signB;
        std::swap(expA, expB);
        std::swap(aSPtr, bSPtr);
    }

    if (0 == expB) {
        expB = 1;

        if (0 == expA) {
            expA = 1;
        }
    }

    uint64_t sigZ = aSPtr->signif;
    uint64_t sigB = bSPtr->signif;
    round_packer_ptr roundPackRoutinePtr = round_pack_to_M_extF80M;
    int32_t const expDiff = expA - expB;

    if (0 != expDiff) {
        uint32_t extSigX[3];
        extSigX[index_word(3, 2)] = sigB >> 32;
        extSigX[index_word(3, 1)] = static_cast<uint32_t>(sigB);
        extSigX[index_word(3, 0)] = 0;
        shift_right_jam_M_96(extSigX, static_cast<uint8_t>(expDiff), extSigX);
        sigB =
            static_cast<uint64_t>(extSigX[index_word(3, 2)]) << 32 |
            extSigX[index_word(3, 1)];

        if (!negateB) {
            sigZ += sigB;

            if (0 == (sigZ & UINT64_C(0x8000000000000000))) {
                auto const sigZExtra_1 = (static_cast<uint32_t>(sigZ) << 31) | !!(0 != extSigX[index_word_lo(3)]);
                ++expA;
                sigZ = UINT64_C(0x8000000000000000) | sigZ >> 1;
                extSigX[index_word(3, 0)] = sigZExtra_1;
            }

            extSigX[index_word(3, 2)] = sigZ >> 32;
            extSigX[index_word(3, 1)] = static_cast<uint32_t>(sigZ);
            (*roundPackRoutinePtr)(signZ, expA, extSigX, extF80_rounding_precision, zSPtr);
            return;
        }

        sigZ -= sigB;
        uint32_t sigZExtra = extSigX[index_word_lo(3)];

        if (0 != sigZExtra) {
            --sigZ;
            sigZExtra = static_cast<uint16_t>(-static_cast<int32_t>(sigZExtra));
        }

        if (0 == (sigZ & UINT64_C(0x8000000000000000))) {
            if (0 != (sigZ & UINT64_C(0x4000000000000000))) {
                --expA;
                sigZ = sigZ << 1 | sigZExtra >> 31;
                sigZExtra <<= 1;
            } else {
                roundPackRoutinePtr = norm_round_pack_to_M_extF80;
            }
        }

        {
            uint32_t extSigX[3];
            extSigX[index_word(3, 0)] = sigZExtra;
            extSigX[index_word(3, 2)] = sigZ >> 32;
            extSigX[index_word(3, 1)] = static_cast<uint32_t>(sigZ);
            (*roundPackRoutinePtr)(signZ, expA, extSigX, extF80_rounding_precision, zSPtr);
        }

        return;
    }

    assert(0 == expDiff);
    uint32_t sigZExtra = 0;

    if (negateB) {
        if (sigZ < sigB) {
            signZ = !signZ;
            sigZ = sigB - sigZ;
        } else {
            sigZ -= sigB;

            if (0 == sigZ) {
                softfloat_round_mode const softfloat_roundingMode = softfloat_get_roundingMode();
                zSPtr->signExp = pack_to_extF80_UI64(softfloat_round_min == softfloat_roundingMode, 0);
                zSPtr->signif = 0;
                return;
            }
        }

        roundPackRoutinePtr = norm_round_pack_to_M_extF80;
    } else {
        sigZ += sigB;

        if (sigZ < sigB) {
            sigZExtra = static_cast<uint32_t>(sigZ) << 31;
            ++expA;
            sigZ = UINT64_C(0x8000000000000000) | sigZ >> 1;
        } else if (0 == (sigZ & UINT64_C(0x8000000000000000))) {
            roundPackRoutinePtr = norm_round_pack_to_M_extF80;
        }
    }

    {
        uint32_t extSigX[3];
        extSigX[index_word(3, 0)] = sigZExtra;
        extSigX[index_word(3, 2)] = sigZ >> 32;
        extSigX[index_word(3, 1)] = static_cast<uint32_t>(sigZ);
        (*roundPackRoutinePtr)(signZ, expA, extSigX, extF80_rounding_precision, zSPtr);
    }
}
#endif

static extFloat80_t
add_magnitudes(extFloat80_t const& a,
               extFloat80_t const& b,
               bool signZ)
{
    using namespace softfloat::internals::fast_int64;

    int32_t const expA = get_exp(a);
    int32_t const expB = get_exp(b);
    int32_t expDiff = expA - expB;

    if (0 == expDiff) {
        if (0x7FFF == expA) {
            extFloat80_t uZ;

            if (0 != (UINT64_C(0x7FFFFFFFFFFFFFFF) & (a.signif | b.signif))) {
                uint128 const uiZ = propagate_NaN(a.signExp, a.signif, b.signExp, b.signif);
                uZ.signExp = static_cast<uint16_t>(uiZ.v64);
                uZ.signif = uiZ.v0;
            } else {
                uZ.signExp = a.signExp;
                uZ.signif = a.signif;
            }

            return uZ;
        }

        auto const sigZ_1 = a.signif + b.signif;

        if (0 == expA) {
            exp32_sig64 const normExpSig = norm_subnormal_extF80Sig(sigZ_1);
            return round_pack_to_extF80(signZ, normExpSig.exp + 1, normExpSig.sig, 0, extF80_rounding_precision);
        }

        uint64_extra const sig64Extra = shortShiftRightJam64Extra(sigZ_1, 0, 1);
        return round_pack_to_extF80(signZ, expA + 1, UINT64_C(0x8000000000000000) | sig64Extra.v, sig64Extra.extra, extF80_rounding_precision);
    }

    if (expDiff < 0) {
        if (0x7FFF == expB) {
            extFloat80_t uZ;

            if (0 != (b.signif & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
                uint128 const uiZ = propagate_NaN(a.signExp, a.signif, b.signExp, b.signif);
                uZ.signExp = static_cast<uint16_t>(uiZ.v64);
                uZ.signif = uiZ.v0;
            } else {
                uZ.signExp = pack_to_extF80_UI64(signZ, 0x7FFF);
                uZ.signif = b.signif;
            }

            return uZ;
        }

        auto const expZ = expB;

        if (0 == expA) {
            ++expDiff;

            if (0 == expDiff) {
                auto const sigZ_1 = a.signif + b.signif;

                if (0 != (UINT64_C(0x8000000000000000) & sigZ_1)) {
                    return round_pack_to_extF80(signZ, expZ, sigZ_1, 0, extF80_rounding_precision);
                }

                uint64_extra const sig64Extra = shortShiftRightJam64Extra(sigZ_1, 0, 1);
                return round_pack_to_extF80(signZ, expZ + 1, UINT64_C(0x8000000000000000) | sig64Extra.v, sig64Extra.extra, extF80_rounding_precision);
            }
        }

        uint64_extra const sig64Extra = shift_right_jam_64Extra(a.signif, 0, static_cast<uint32_t>(-expDiff));
        auto const sigZ_1 = sig64Extra.v + b.signif;

        if (0 != (sigZ_1 & UINT64_C(0x8000000000000000))) {
            return round_pack_to_extF80(signZ, expZ, sigZ_1, sig64Extra.extra, extF80_rounding_precision);
        }

        auto const sig64Extra_1 = shortShiftRightJam64Extra(sigZ_1, sig64Extra.extra, 1);
        return round_pack_to_extF80(signZ, expZ + 1, sig64Extra_1.v | UINT64_C(0x8000000000000000), sig64Extra_1.extra, extF80_rounding_precision);
    }

    if (0x7FFF == expA) {
        if (0 != (UINT64_C(0x7FFFFFFFFFFFFFFF) & a.signif)) {
            uint128 const uiZ = propagate_NaN(a.signExp, a.signif, b.signExp, b.signif);
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(uiZ.v64);
            uZ.signif = uiZ.v0;
            return uZ;
        } else {
            extFloat80_t uZ;
            uZ.signExp = a.signExp;
            uZ.signif = a.signif;
            return uZ;
        }
    }

    auto const expZ = expA;

    if (0 == expB) {
        --expDiff;

        if (0 == expDiff) {
            auto const sigZ_1 = a.signif + b.signif;

            if (0 != (UINT64_C(0x8000000000000000) & sigZ_1)) {
                return round_pack_to_extF80(signZ, expZ, sigZ_1, 0, extF80_rounding_precision);
            }

            uint64_extra const sig64Extra = shortShiftRightJam64Extra(sigZ_1, 0, 1);
            return round_pack_to_extF80(signZ, expZ + 1, UINT64_C(0x8000000000000000) | sig64Extra.v, sig64Extra.extra, extF80_rounding_precision);
        }
    }

    uint64_extra const sig64Extra = shift_right_jam_64Extra(b.signif, 0u, static_cast<uint32_t>(expDiff));
    auto const sigZ_1 = a.signif + sig64Extra.v;

    if (0 != (sigZ_1 & UINT64_C(0x8000000000000000))) {
        return round_pack_to_extF80(signZ, expZ, sigZ_1, sig64Extra.extra, extF80_rounding_precision);
    }

    auto const sig64Extra_1 = shortShiftRightJam64Extra(sigZ_1, sig64Extra.extra, 1);
    return round_pack_to_extF80(signZ, expZ + 1, sig64Extra_1.v | UINT64_C(0x8000000000000000), sig64Extra_1.extra, extF80_rounding_precision);
}

static extFloat80_t
sub_magnitudes(extFloat80_t const& a,
               extFloat80_t const& b,
               bool const signZ)
{
    using namespace softfloat::internals::fast_int64;

    int32_t const expA = exp_extF80_UI64(a.signExp);
    int32_t const expB = exp_extF80_UI64(b.signExp);
    int32_t expDiff = expA - expB;

    if (0 < expDiff) {
        if (0x7FFF == expA) {
            if (0 != (a.signif & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
                uint128 const uiZ = propagate_NaN(a.signExp, a.signif, b.signExp, b.signif);
                extFloat80_t uZ;
                uZ.signExp = static_cast<uint16_t>(uiZ.v64);
                uZ.signif = uiZ.v0;
                return uZ;
            }

            extFloat80_t uZ;
            uZ.signExp = a.signExp;
            uZ.signif = a.signif;
            return uZ;
        }

        if (0 == expB) {
            --expDiff;

            if (0 == expDiff) {
                uint128 const sig128_6 = sub(uint128{a.signif, 0}, uint128{b.signif, 0});
                return
                    norm_round_pack_to_extF80(signZ, expA, sig128_6.v64, sig128_6.v0, extF80_rounding_precision);
            }
        }

        uint128 const sig128_2 = shift_right_jam_128(b.signif, 0, static_cast<uint32_t>(expDiff));
        uint128 const sig128_1 = sub(uint128{a.signif, 0}, sig128_2);
        return
            norm_round_pack_to_extF80(signZ, expA, sig128_1.v64, sig128_1.v0, extF80_rounding_precision);
    }

    if (expDiff < 0) {
        if (expB == 0x7FFF) {
            if (0 != (b.signif & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
                uint128 const uiZ = propagate_NaN(a.signExp, a.signif, b.signExp, b.signif);
                extFloat80_t uZ;
                uZ.signExp = static_cast<uint16_t>(uiZ.v64);
                uZ.signif = uiZ.v0;
                return uZ;
            }

            extFloat80_t uZ;
            uZ.signExp = pack_to_extF80_UI64(!signZ, 0x7FFF);
            uZ.signif = UINT64_C(0x8000000000000000);
            return uZ;
        }

        if (0 == expA) {
            ++expDiff;

            if (0 == expDiff) {
                uint128 const sig128_7 = sub(uint128{b.signif, 0}, uint128{a.signif, 0});
                return
                    norm_round_pack_to_extF80(!signZ, expB, sig128_7.v64, sig128_7.v0, extF80_rounding_precision);
            }
        }

        uint128 const sig128_4 = shift_right_jam_128(a.signif, 0, static_cast<uint32_t>(-expDiff));
        uint128 const sig128_3 = sub(uint128{b.signif, 0}, sig128_4);
        return
            norm_round_pack_to_extF80(!signZ, expB, sig128_3.v64, sig128_3.v0, extF80_rounding_precision);
    }

    if (0x7FFF == expA) {
        if (0 != ((a.signif | b.signif) & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
            uint128 const uiZ = propagate_NaN(a.signExp, a.signif, b.signExp, b.signif);
            extFloat80_t uZ;
            uZ.signExp = static_cast<uint16_t>(uiZ.v64);
            uZ.signif = uiZ.v0;
            return uZ;
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        extFloat80_t uZ;
        uZ.signExp = defaultNaNExtF80UI64;
        uZ.signif = defaultNaNExtF80UI0;
        return uZ;
    }

    int32_t const expZ = 0 == expA ? 1 : expA;

    if (a.signif > b.signif) {
        uint128 const sig128 = sub(uint128{a.signif, 0}, uint128{b.signif, 0});
        return norm_round_pack_to_extF80(signZ, expZ, sig128.v64, sig128.v0, extF80_rounding_precision);
    }

    if (a.signif < b.signif) {
        uint128 const sig128 = sub(uint128{b.signif, 0}, uint128{a.signif, 0});
        return norm_round_pack_to_extF80(!signZ, expZ, sig128.v64, sig128.v0, extF80_rounding_precision);
    }

    /* a.signif == b.signif */
    extFloat80_t uZ;
    uZ.signExp = pack_to_extF80_UI64(softfloat_round_min == softfloat_get_roundingMode(), 0);
    uZ.signif = 0;
    return uZ;
}

}  // namespace

extFloat80_t
extF80_add(extFloat80_t const a,
           extFloat80_t const b)
{
    using namespace softfloat::internals::fast_int64;
    bool const signA = is_sign(a);
    bool const signB = is_sign(b);
    return signA == signB ?
        add_magnitudes(a, b, signA) :
        sub_magnitudes(a, b, signA);
}

extFloat80_t
extF80_sub(extFloat80_t const a,
           extFloat80_t const b)
{
    using namespace softfloat::internals::fast_int64;

    bool const signA = is_sign(a);
    bool const signB = is_sign(b);
    return
        signA == signB ?
        sub_magnitudes(a, b, signA) :
        add_magnitudes(a, b, signA);
}

void
extF80M_add(extFloat80_t const *const aPtr,
            extFloat80_t const *const bPtr,
            extFloat80_t *const zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    using namespace softfloat::internals::fast_int64;

    bool const signA = is_sign(*aPtr);
    bool const signB = is_sign(*bPtr);
    *zPtr =
        signA == signB ?
        add_magnitudes(*aPtr, *bPtr, signA) :
        sub_magnitudes(*aPtr, *bPtr, signA);
#else

    using namespace softfloat::internals::slow_int64;
    add_M_extF80(aPtr, bPtr, zPtr, false);

#endif
}

void
extF80M_sub(extFloat80_t const *const aPtr,
            extFloat80_t const *const bPtr,
            extFloat80_t* const zPtr)
{
#if (SOFTFLOAT_FAST_INT64)

    using namespace softfloat::internals::fast_int64;

    bool const signA = is_sign(*aPtr);
    bool const signB = is_sign(*bPtr);
    *zPtr =
        signA == signB ?
        sub_magnitudes(*aPtr, *bPtr, signA) :
        add_magnitudes(*aPtr, *bPtr, signA);

#else

    using namespace softfloat::internals::slow_int64;
    add_M_extF80(aPtr, bPtr, zPtr, true);

#endif
}
