
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

#if (SOFTFLOAT_FAST_INT64)

void
f128M_rem(float128_t const* const aPtr,
          float128_t const* const bPtr,
          float128_t* const zPtr)
{
    *zPtr = f128_rem(*aPtr, *bPtr);
}

#else

using namespace softfloat::internals::slow_int64;

namespace {
static inline void
copyA(uint32_t const uiA96,
      uint32_t const* const aWPtr,
      uint32_t zWPtr[4])
{
    zWPtr[index_word_hi(4)] = uiA96;
    zWPtr[index_word(4, 2)] = aWPtr[index_word(4, 2)];
    zWPtr[index_word(4, 1)] = aWPtr[index_word(4, 1)];
    zWPtr[index_word(4, 0)] = aWPtr[index_word(4, 0)];
}
}  // namespace

void
f128M_rem(float128_t const* const aPtr,
          float128_t const* const bPtr,
          float128_t* const zPtr)
{
    auto const aWPtr = reinterpret_cast<const uint32_t*>(aPtr);
    auto const bWPtr = reinterpret_cast<const uint32_t*>(bPtr);
    auto const zWPtr = reinterpret_cast<uint32_t*>(zPtr);

    uint32_t const uiA96 = aWPtr[index_word_hi(4)];
    int32_t expA = exp_F128_UI96(uiA96);
    int32_t expB = exp_F128_UI96(bWPtr[index_word_hi(4)]);

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (try_propagate_NaN_M_F128(aWPtr, bWPtr, zWPtr)) {
            return;
        }

        if (expA == 0x7FFF) {
            invalid_M_F128(zWPtr);
            return;
        }

        copyA(uiA96, aWPtr, zWPtr);
        return;
    }

    if (expA < expB - 1) {
        copyA(uiA96, aWPtr, zWPtr);
        return;
    }

    uint32_t x[4];
    expB = shift_norm_sig_M_F128(bWPtr, 13, x);

    if (-128 == expB) {
        invalid_M_F128(zWPtr);
        return;
    }

    uint32_t rem1[5];
    uint32_t* remPtr = &rem1[index_multiword_lo(5, 4)];
    expA = shift_norm_sig_M_F128(aWPtr, 13, remPtr);

    if (expA == -128) {
        copyA(uiA96, aWPtr, zWPtr);
        return;
    }

    bool signRem = is_sign(uiA96);

    int32_t expDiff = expA - expB;
    uint32_t* altRemPtr;
    uint32_t rem2[5];
    uint32_t q;

    if (expDiff < 1) {
        if (expDiff < -1) {
            copyA(uiA96, aWPtr, zWPtr);
            return;
        }

        if (expDiff) {
            --expB;
            add_M_128(x, x, x);
            q = 0;
        } else {
            q = 0u + !!(compare_M_128(x, remPtr) <= 0);

            if (q) {
                sub_M_128(remPtr, x, remPtr);
            }
        }
    } else {
        uint32_t const recip32 =
            approx_recip_32_1((static_cast<uint64_t>(x[index_word(4, 3)]) << 32 | x[index_word(4, 2)]) >> 30);
        expDiff -= 30;

        uint64_t q64;

        for (;;) {
            q64 = static_cast<uint64_t>(remPtr[index_word_hi(4)]) * recip32;

            if (expDiff < 0) {
                break;
            }

            q = (q64 + 0x80000000) >> 32;
            rem_step_by_32_M_128(remPtr, 29, x, q, remPtr);

            if (remPtr[index_word_hi(4)] & 0x80000000) {
                add_M_128(remPtr, x, remPtr);
            }

            expDiff -= 29;
        }

        /* `expDiff' cannot be less than -29 here. */
        q = static_cast<uint32_t>(q64 >> 32) >> (~expDiff & 31);
        /**
        @todo Warning   C4244   '=': conversion from 'int' to 'uint8_t', possible loss of data
        */
        rem_step_by_32_M_128(remPtr, static_cast<uint8_t>(expDiff + 30), x, q, remPtr);

        if (0 != (remPtr[index_word_hi(4)] & 0x80000000)) {
            altRemPtr = &rem2[index_multiword_lo(5, 4)];
            add_M_128(remPtr, x, altRemPtr);
            goto selectRem;
        }
    }

    altRemPtr = &rem2[index_multiword_lo(5, 4)];

    do {
        ++q;
        uint32_t* const newRemPtr = altRemPtr;
        sub_M_128(remPtr, x, newRemPtr);
        altRemPtr = remPtr;
        remPtr = newRemPtr;
    } while (0 == (remPtr[index_word_hi(4)] & 0x80000000));

selectRem:

    {
        add_M_128(remPtr, altRemPtr, x);

        if (
            0 != (x[index_word_hi(4)] & 0x80000000) ||
            (
                0 != (q & 1) &&
                0 == x[index_word(4, 0)] &&
                0 == x[index_word(4, 1)] &&
                0 == x[index_word(4, 2)] &&
                0 == x[index_word_hi(4)]
            )
        ) {
            remPtr = altRemPtr;
        }

        if (0 != (remPtr[index_word_hi(4)] & 0x80000000)) {
            signRem = !signRem;
            neg_M_X128(remPtr);
        }

        remPtr -= index_multiword_lo(5, 4);
        remPtr[index_word_hi(5)] = 0;
        norm_round_pack_to_M_F128(signRem, expB + 18, remPtr, zWPtr);
        return;
    }
}
#endif
