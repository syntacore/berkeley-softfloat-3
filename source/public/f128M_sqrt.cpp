
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
f128M_sqrt(float128_t const* const aPtr,
           float128_t* const zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    *zPtr = f128_sqrt(*aPtr);
#else
    using namespace softfloat::internals::slow_int64;

    auto const aWPtr = reinterpret_cast<uint32_t const*>(aPtr);
    auto const zWPtr = reinterpret_cast<uint32_t*>(zPtr);
    uint32_t const uiA96 = aWPtr[index_word_hi(4)];
    bool const signA = is_sign(uiA96);
    int32_t const rawExpA = exp_F128_UI96(uiA96);

    if (0x7FFF == rawExpA) {
        if (frac_F128_UI96(uiA96) || 0 != (aWPtr[index_word(4, 2)] | aWPtr[index_word(4, 1)] | aWPtr[index_word(4, 0)])) {
            propagate_NaN_F128M(aWPtr, 0, zWPtr);
            return;
        }

        if (!signA) {
            zWPtr[index_word_hi(4)] = uiA96;
            zWPtr[index_word(4, 2)] = aWPtr[index_word(4, 2)];
            zWPtr[index_word(4, 1)] = aWPtr[index_word(4, 1)];
            zWPtr[index_word(4, 0)] = aWPtr[index_word(4, 0)];
            return;
        }

        invalid_M_F128(zWPtr);
        return;
    }

    uint32_t rem[6];
    int32_t expA = shift_norm_sig_M_F128(aWPtr, static_cast<uint8_t>(13 - (rawExpA & 1)), rem);

    if (-128 == expA) {
        zWPtr[index_word_hi(4)] = uiA96;
        zWPtr[index_word(4, 2)] = aWPtr[index_word(4, 2)];
        zWPtr[index_word(4, 1)] = aWPtr[index_word(4, 1)];
        zWPtr[index_word(4, 0)] = aWPtr[index_word(4, 0)];
        return;
    }

    if (signA) {
        invalid_M_F128(zWPtr);
        return;
    }

    /*
    `sig32Z' is guaranteed to be a lower bound on the square root of
    `sig32A', which makes `sig32Z' also a lower bound on the square root of
    `sigA'.
    */
    int32_t const expZ = ((expA - 0x3FFF) >> 1) + 0x3FFE;
    expA &= 1;
    uint64_t rem64 = static_cast<uint64_t>(rem[index_word(4, 3)]) << 32 | rem[index_word(4, 2)];
    uint32_t sig32A;

    if (expA) {
        if (!rawExpA) {
            short_shift_right_M_128(rem, 1, rem);
            rem64 >>= 1;
        }

        sig32A = static_cast<uint32_t>(rem64 >> 29);
    } else {
        sig32A = static_cast<uint32_t>(rem64 >> 30);
    }

    uint32_t const recipSqrt32 = approx_recip_sqrt_32_1(static_cast<uint32_t>(expA), sig32A);
    uint32_t sig32Z = (static_cast<uint64_t>(sig32A) * recipSqrt32) >> 32;

    if (0 != expA) {
        sig32Z >>= 1;
    }

    uint32_t qs[3];
    qs[2] = sig32Z;
    rem64 -= static_cast<uint64_t>(sig32Z) * sig32Z;
    rem[index_word(4, 3)] = rem64 >> 32;
    rem[index_word(4, 2)] = static_cast<uint32_t>(rem64);

    uint32_t q = (static_cast<uint32_t>(rem64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    qs[1] = q;
    uint64_t sig64Z = (static_cast<uint64_t>(sig32Z) << 32) + (static_cast<uint64_t>(q) << 3);
    uint64_t x64 = (static_cast<uint64_t>(sig32Z) << 32) + sig64Z;
    uint32_t term[5];
    term[index_word(4, 3)] = 0;
    term[index_word(4, 2)] = x64 >> 32;
    term[index_word(4, 1)] = static_cast<uint32_t>(x64);
    term[index_word(4, 0)] = 0;
    uint32_t y[5];
    rem_step_by_32_M_128(rem, 29, term, q, y);
    rem64 = static_cast<uint64_t>(y[index_word(4, 3)]) << 32 | y[index_word(4, 2)];

    q = (static_cast<uint32_t>(rem64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    sig64Z <<= 1;

    /* Repeating this loop is a rare occurrence. */
    uint32_t rem32;

    for (;;) {
        x64 = sig64Z + (q >> 26);
        term[index_word(4, 2)] = x64 >> 32;
        term[index_word(4, 1)] = static_cast<uint32_t>(x64);
        term[index_word(4, 0)] = q << 6;
        term[index_word(4, 3)] = 0;
        rem_step_by_32_M_128(y, 29, term, q, &rem[index_multiword_hi(6, 4)]);
        rem32 = rem[index_word_hi(6)];

        if (!(rem32 & 0x80000000)) {
            break;
        }

        --q;
    }

    qs[0] = q;
    rem64 = static_cast<uint64_t>(rem32) << 32 | rem[index_word(6, 4)];

    q = ((static_cast<uint32_t>(rem64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32) + 2;
    x64 = static_cast<uint64_t>(q) << 27;
    y[index_word(5, 0)] = static_cast<uint32_t>(x64);
    x64 = (static_cast<uint64_t>(qs[0]) << 24) + (x64 >> 32);
    y[index_word(5, 1)] = static_cast<uint32_t>(x64);
    x64 = (static_cast<uint64_t>(qs[1]) << 21) + (x64 >> 32);
    y[index_word(5, 2)] = static_cast<uint32_t>(x64);
    x64 = (static_cast<uint64_t>(qs[2]) << 18) + (x64 >> 32);
    y[index_word(5, 3)] = static_cast<uint32_t>(x64);
    y[index_word(5, 4)] = x64 >> 32;

    if ((q & 0xF) <= 2) {
        q &= ~3;
        y[index_word_lo(5)] = q << 27;
        term[index_word(5, 4)] = 0;
        term[index_word(5, 3)] = 0;
        term[index_word(5, 2)] = 0;
        term[index_word(5, 1)] = q >> 6;
        term[index_word(5, 0)] = q << 26;
        sub_M_160(y, term, term);
        rem[index_word(6, 1)] = 0;
        rem[index_word(6, 0)] = 0;
        rem_step_M_160_by_32(
            &rem[index_multiword_lo(6, 5)],
            14,
            term,
            q,
            &rem[index_multiword_lo(6, 5)]
        );
        rem32 = rem[index_word(6, 4)];

        if (rem32 & 0x80000000) {
            sub_1_M_X160(y);
        } else if (0 != rem32 || 0 != rem[index_word(6, 0)] || 0 != rem[index_word(6, 1)] || 0 != (rem[index_word(6, 3)] | rem[index_word(6, 2)])) {
            y[index_word_lo(5)] |= 1;
        }
    }

    round_pack_to_M_F128(0, expZ, y, zWPtr);
    return;
#endif
}
