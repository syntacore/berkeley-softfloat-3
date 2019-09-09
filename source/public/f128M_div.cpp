
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

void
f128M_div(float128_t const* const aPtr,
          float128_t const* const bPtr,
          float128_t* const zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    *zPtr = f128_div(*aPtr, *bPtr);
#else
    using namespace softfloat::internals::slow_int64;

    auto const aWPtr = reinterpret_cast<uint32_t const*>(aPtr);
    auto const bWPtr = reinterpret_cast<uint32_t const*>(bPtr);
    auto const zWPtr = reinterpret_cast<uint32_t*>(zPtr);

    uint32_t const uiA96 = aWPtr[index_word_hi(4)];
    bool const signA = is_sign(uiA96);
    int32_t expA = exp_F128_UI96(uiA96);
    uint32_t const uiB96 = bWPtr[index_word_hi(4)];
    bool const signB = is_sign(uiB96);
    int32_t expB = exp_F128_UI96(uiB96);
    bool const signZ = signA != signB;

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (try_propagate_NaN_M_F128(aWPtr, bWPtr, zWPtr)) {
            return;
        }

        if (expA == 0x7FFF) {
            if (expB == 0x7FFF) {
                invalid_M_F128(zWPtr);
                return;
            }

            zWPtr[index_word_hi(4)] = pack_to_F128_UI96(signZ, 0x7FFF, 0);
            zWPtr[index_word(4, 2)] = 0;
            zWPtr[index_word(4, 1)] = 0;
            zWPtr[index_word(4, 0)] = 0;
            return;
        }

        zWPtr[index_word_hi(4)] = pack_to_F128_UI96(signZ, 0, 0);
        zWPtr[index_word(4, 2)] = 0;
        zWPtr[index_word(4, 1)] = 0;
        zWPtr[index_word(4, 0)] = 0;
        return;
    }

    uint32_t y[5];
    expA = shift_norm_sig_M_F128(aWPtr, 13, y);
    uint32_t sigB[4];
    expB = shift_norm_sig_M_F128(bWPtr, 13, sigB);

    if (expA == -128) {
        if (expB == -128) {
            invalid_M_F128(zWPtr);
            return;
        }

        zWPtr[index_word_hi(4)] = pack_to_F128_UI96(signZ, 0, 0);
        zWPtr[index_word(4, 2)] = 0;
        zWPtr[index_word(4, 1)] = 0;
        zWPtr[index_word(4, 0)] = 0;
        return;
    }

    if (expB == -128) {
        softfloat_raiseFlags(softfloat_flag_infinite);
        zWPtr[index_word_hi(4)] = pack_to_F128_UI96(signZ, 0x7FFF, 0);
        zWPtr[index_word(4, 2)] = 0;
        zWPtr[index_word(4, 1)] = 0;
        zWPtr[index_word(4, 0)] = 0;
        return;
    }

    int32_t expZ = expA - expB + 0x3FFE;

    if (compare_M_128(y, sigB) < 0) {
        --expZ;
        add_M_128(y, y, y);
    }

    uint32_t const recip32 = approx_recip_32_1((static_cast<uint64_t>(sigB[index_word(4, 3)]) << 32 | sigB[index_word(4, 2)]) >> 30);
    int ix = 3;
    uint64_t q64;
    uint32_t qs[3];
    uint32_t q;

    for (;;) {
        q64 = static_cast<uint64_t>(y[index_word_hi(4)]) * recip32;
        q = (q64 + 0x80000000) >> 32;
        --ix;

        if (ix < 0) {
            break;
        }

        rem_step_by_32_M_128(y, 29, sigB, q, y);

        if (y[index_word_hi(4)] & 0x80000000) {
            --q;
            add_M_128(y, sigB, y);
        }

        qs[ix] = q;
    }

    if (((q + 1) & 7) < 2) {
        rem_step_by_32_M_128(y, 29, sigB, q, y);

        if (y[index_word_hi(4)] & 0x80000000) {
            --q;
            add_M_128(y, sigB, y);
        } else if (compare_M_128(sigB, y) <= 0) {
            ++q;
            sub_M_128(y, sigB, y);
        }

        if (y[index_word_lo(4)] || y[index_word(4, 1)] || (y[index_word(4, 2)] | y[index_word(4, 3)])) {
            q |= 1;
        }
    }

    q64 = static_cast<uint64_t>(q) << 28;
    y[index_word(5, 0)] = static_cast<uint32_t>(q64);
    q64 = (static_cast<uint64_t>(qs[0]) << 25) + (q64 >> 32);
    y[index_word(5, 1)] = static_cast<uint32_t>(q64);
    q64 = (static_cast<uint64_t>(qs[1]) << 22) + (q64 >> 32);
    y[index_word(5, 2)] = static_cast<uint32_t>(q64);
    q64 = (static_cast<uint64_t>(qs[2]) << 19) + (q64 >> 32);
    y[index_word(5, 3)] = static_cast<uint32_t>(q64);
    y[index_word(5, 4)] = q64 >> 32;
    round_pack_to_M_F128(signZ, expZ, y, zWPtr);
    return;
#endif
}
