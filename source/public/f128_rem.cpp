
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

float128_t
f128_rem(float128_t const a,
         float128_t const b)
{
    using namespace softfloat::internals::fast_int64;

    uint128 const uA{a};
    bool const signA = is_sign(uA.v64);
    int32_t expA = exp_F128_UI64(uA.v64);
    uint128 sigA{frac_F128_UI64(uA.v64), uA.v0};
    uint128 const uB{b};
    int32_t expB = exp_F128_UI64(uB.v64);
    uint128 sigB{frac_F128_UI64(uB.v64), uB.v0};

    if (0x7FFF == expA) {
        if (0 != (sigA.v64 | sigA.v0) || (0x7FFF == expB && 0 != (sigB.v64 | sigB.v0))) {
            return float128_t(propagate_NaN(uA, uB));
        }

        softfloat_raiseFlags(softfloat_flag_invalid);
        return static_cast<float128_t>(uint128{defaultNaNF128UI64, defaultNaNF128UI0});
    }

    if (0x7FFF == expB) {
        return
            0 == (sigB.v64 | sigB.v0) ?
            a :
            float128_t(propagate_NaN(uA, uB));
    }

    if (0 == expB) {
        if (0 == (sigB.v64 | sigB.v0)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return static_cast<float128_t>(uint128{defaultNaNF128UI64, defaultNaNF128UI0});
        }

        exp32_sig128 const normExpSig = norm_subnormal_F128Sig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    if (0 == expA) {
        if (0 == (sigA.v64 | sigA.v0)) {
            return a;
        }

        exp32_sig128 const normExpSig = norm_subnormal_F128Sig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    sigA.v64 |= UINT64_C(0x0001000000000000);
    sigB.v64 |= UINT64_C(0x0001000000000000);

    uint128 rem = sigA;
    int32_t expDiff = expA - expB;
    uint32_t q;
    uint128 altRem;

    if (expDiff < 1) {
        if (expDiff < -1) {
            return a;
        }

        if (expDiff) {
            --expB;
            sigB = add(sigB, sigB);
            q = 0;
        } else {
            /* suppress warning: '=': conversion from 'bool' to 'uint32_t', signed/unsigned mismatch */
            q = 0u + !!le(sigB, rem);

            if (0 != q) {
                rem = sub(rem, sigB);
            }
        }
    } else {
        uint32_t const recip32 = approx_recip_32_1(static_cast<uint32_t>(sigB.v64 >> 17));
        expDiff -= 30;
        uint64_t q64;

        for (;;) {
            q64 = static_cast<uint64_t>(static_cast<uint32_t>(rem.v64 >> 19)) * recip32;

            if (expDiff < 0) {
                break;
            }

            q = (q64 + 0x80000000) >> 32;
            rem = short_shift_left_128(rem, 29);
            rem = sub(rem, mul_128_by_32(sigB, q));

            if (0 != (rem.v64 & UINT64_C(0x8000000000000000))) {
                rem = add(rem, sigB);
            }

            expDiff -= 29;
        }

        /* `expDiff' cannot be less than -29 here.*/
        assert(-29 <= expDiff);
        q = static_cast<uint32_t>(q64 >> 32) >> (~expDiff & 31);
        rem = short_shift_left_128(rem, static_cast<uint8_t>(expDiff + 30u));
        rem = sub(rem, mul_128_by_32(sigB, q));

        if (0 != (rem.v64 & UINT64_C(0x8000000000000000))) {
            altRem = add(rem, sigB);
            uint128 const meanRem = add(rem, altRem);

            if (0 != (meanRem.v64 & UINT64_C(0x8000000000000000)) || (0 == (meanRem.v64 | meanRem.v0) && 0 != (q & 1))) {
                rem = altRem;
            }

            if (0 != (rem.v64 & UINT64_C(0x8000000000000000))) {
                auto const rem_1 = sub(uint128{0, 0}, rem);
                return norm_round_pack_to_F128(!signA, expB - 1, rem_1.v64, rem_1.v0);
            }

            return norm_round_pack_to_F128(signA, expB - 1, rem.v64, rem.v0);
        }
    }

    do {
        altRem = rem;
        ++q;
        rem = sub(rem, sigB);
    } while (0 == (rem.v64 & UINT64_C(0x8000000000000000)));

    uint128 const meanRem = add(rem, altRem);

    if (0 != (meanRem.v64 & UINT64_C(0x8000000000000000)) || (0 == (meanRem.v64 | meanRem.v0) && 0 != (q & 1))) {
        rem = altRem;
    }

    if (0 != (rem.v64 & UINT64_C(0x8000000000000000))) {
        auto const rem_1 = sub(uint128{0, 0}, rem);
        return norm_round_pack_to_F128(!signA, expB - 1, rem_1.v64, rem_1.v0);
    }

    return norm_round_pack_to_F128(signA, expB - 1, rem.v64, rem.v0);
}
