
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All rights reserved.
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

void
extF80M_div(extFloat80_t const* const aSPtr,
            extFloat80_t const* const bSPtr,
            extFloat80_t* const zSPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    *zSPtr = extF80_div(*aSPtr, *bSPtr);
#else
    using namespace softfloat::internals::slow_int64;

    int32_t expA = get_exp(*aSPtr);
    int32_t expB = get_exp(*bSPtr);
    bool const signZ = is_sign(*aSPtr) != is_sign(*bSPtr);

    if (expA == INT16_MAX || expB == INT16_MAX) {
        if (try_propagate_NaN_M_extF80(aSPtr, bSPtr, zSPtr)) {
            return;
        }

        if (expA != INT16_MAX) {
            zSPtr->signExp = pack_to_extF80_UI64(signZ, 0);
            zSPtr->signif = 0;
            return;
        }

        if (expB == INT16_MAX) {
            invalid_M_extF80(zSPtr);
            return;
        }

        zSPtr->signExp = pack_to_extF80_UI64(signZ, 0x7FFF);
        zSPtr->signif = UINT64_C(0x8000000000000000);
        return;
    }

    uint64_t sigA = aSPtr->signif;
    uint64_t x64 = bSPtr->signif;

    if (!expB) {
        expB = 1;
    }

    if (0 == (x64 & UINT64_C(0x8000000000000000))) {
        if (!x64) {
            if (0 == sigA) {
                invalid_M_extF80(zSPtr);
                return;
            }

            softfloat_raiseFlags(softfloat_flag_infinite);
            zSPtr->signExp = pack_to_extF80_UI64(signZ, 0x7FFF);
            zSPtr->signif = UINT64_C(0x8000000000000000);
            return;
        }

        expB += norm_M_extF80Sig(&x64);
    }

    if (!expA) {
        expA = 1;
    }

    if (0 == (sigA & UINT64_C(0x8000000000000000))) {
        if (0 == sigA) {
            zSPtr->signExp = pack_to_extF80_UI64(signZ, 0);
            zSPtr->signif = 0;
            return;
        }

        expA += norm_M_extF80Sig(&sigA);
    }

    int32_t expZ = expA - expB + 0x3FFF;
    uint8_t shiftDist = 29;

    if (sigA < x64) {
        --expZ;
        shiftDist = 30;
    }

    uint32_t y[3];
    short_shift_left_M_64_to_96(sigA, shiftDist, y);
    uint32_t const recip32 = approx_recip_32_1(static_cast<uint32_t>(x64 >> 32));
    uint32_t sigB[3];
    sigB[index_word(3, 0)] = static_cast<uint32_t>(x64) << 30;
    x64 >>= 2;
    sigB[index_word(3, 2)] = x64 >> 32;
    sigB[index_word(3, 1)] = static_cast<uint32_t>(x64);
    int ix = 2;
    uint32_t q;
    uint32_t qs[2];

    for (;;) {
        x64 = static_cast<uint64_t>(y[index_word_hi(3)]) * recip32;
        q = (x64 + 0x80000000) >> 32;
        --ix;

        if (ix < 0) {
            break;
        }

        rem_step_by_32_M_96(y, 29, sigB, q, y);

        if (y[index_word_hi(3)] & 0x80000000) {
            --q;
            add_M_96(y, sigB, y);
        }

        qs[ix] = q;
    }

    if (((q + 1) & 0x3FFFFF) < 2) {
        rem_step_by_32_M_96(y, 29, sigB, q, y);

        if (y[index_word_hi(3)] & 0x80000000) {
            --q;
            add_M_96(y, sigB, y);
        } else if (compare_M_96(sigB, y) <= 0) {
            ++q;
            sub_M_96(y, sigB, y);
        }

        if (y[index_word_lo(3)] || y[index_word(3, 1)] || y[index_word(3, 2)]) {
            q |= 1;
        }
    }

    x64 = static_cast<uint64_t>(q) << 9;
    y[index_word(3, 0)] = static_cast<uint32_t>(x64);
    x64 = (static_cast<uint64_t>(qs[0]) << 6) + (x64 >> 32);
    y[index_word(3, 1)] = static_cast<uint32_t>(x64);
    y[index_word(3, 2)] = (qs[1] << 3) + static_cast<uint32_t>(x64 >> 32);
    round_pack_to_M_extF80M(signZ, expZ, y, extF80_rounding_precision, zSPtr);
    return;
#endif
}
