
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015 The Regents of the University of
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
extF80M_sqrt(extFloat80_t const* const aPtr, extFloat80_t* const zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    *zPtr = extF80_sqrt(*aPtr);
#else
    using namespace softfloat::internals::slow_int64;
    uint16_t uiA64 = aPtr->signExp;
    uint16_t const signUI64 = static_cast<uint16_t>(uiA64 & pack_to_extF80_UI64(true, 0u));
    int32_t expA = exp_extF80_UI64(uiA64);
    uint64_t rem64 = aPtr->signif;

    if (expA == 0x7FFF) {
        if (rem64 & UINT64_C(0x7FFFFFFFFFFFFFFF)) {
            propagate_NaN_extF80M(aPtr, 0, zPtr);
            return;
        }

        if (signUI64) {
            invalid_M_extF80(zPtr);
            return;
        }

        rem64 = UINT64_C(0x8000000000000000);
        zPtr->signExp = uiA64;
        zPtr->signif = rem64;
        return;
    }

    if (!expA) {
        expA = 1;
    }

    if (!(rem64 & UINT64_C(0x8000000000000000))) {
        if (!rem64) {
            zPtr->signExp = signUI64;
            zPtr->signif = rem64;
            return;
        }

        expA += norm_M_extF80Sig(&rem64);
    }

    if (signUI64) {
        invalid_M_extF80(zPtr);
        return;
    }

    int32_t const expZ = ((expA - 0x3FFF) >> 1) + 0x3FFF;
    expA &= 1;
    uint32_t rem[4];
    short_shift_left_M_64_to_96(rem64, static_cast<uint8_t>(30 - expA), &rem[index_multiword_hi(4, 3)]);
    uint32_t const sig32A = rem64 >> 32;
    uint32_t const recipSqrt32 = approx_recip_sqrt_32_1(static_cast<uint32_t>(expA), sig32A);
    uint32_t sig32Z = (static_cast<uint64_t>(sig32A) * recipSqrt32) >> 32;

    if (expA) {
        sig32Z >>= 1;
    }

    rem64 = 
        (static_cast<uint64_t>(rem[index_word(4, 3)]) << 32 | rem[index_word(4, 2)]) - 
        static_cast<uint64_t>(sig32Z) * sig32Z;
    rem[index_word(4, 3)] = rem64 >> 32;
    rem[index_word(4, 2)] = static_cast<uint32_t>(rem64);

    uint32_t q = (static_cast<uint32_t>(rem64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32;
    uint64_t const sig64Z = (static_cast<uint64_t>(sig32Z) << 32) + (static_cast<uint64_t>(q) << 3);
    uint64_t x64 = (static_cast<uint64_t>(sig32Z) << 32) + sig64Z;
    uint32_t term[4];
    term[index_word(3, 2)] = 0;
    term[index_word(3, 1)] = x64 >> 32;
    term[index_word(3, 0)] = static_cast<uint32_t>(x64);
    rem_step_by_32_M_96(
        &rem[index_multiword_hi(4, 3)],
        29,
        term,
        q,
        &rem[index_multiword_hi(4, 3)]
    );
    rem64 = static_cast<uint64_t>(rem[index_word(4, 3)]) << 32 | rem[index_word(4, 2)];

    q = ((static_cast<uint32_t>(rem64 >> 2) * static_cast<uint64_t>(recipSqrt32)) >> 32) + 2;
    x64 = static_cast<uint64_t>(q) << 7;
    uint32_t extSigZ[3];
    extSigZ[index_word(3, 0)] = static_cast<uint32_t>(x64);
    x64 = (sig64Z << 1) + (x64 >> 32);
    extSigZ[index_word(3, 2)] = x64 >> 32;
    extSigZ[index_word(3, 1)] = static_cast<uint32_t>(x64);

    if ((q & 0xFFFFFF) <= 2) {
        q &= ~static_cast<uint32_t>(0xFFFF);
        extSigZ[index_word_lo(3)] = q << 7;
        x64 = sig64Z + (q >> 27);
        term[index_word(4, 3)] = 0;
        term[index_word(4, 2)] = x64 >> 32;
        term[index_word(4, 1)] = static_cast<uint32_t>(x64);
        term[index_word(4, 0)] = q << 5;
        rem[index_word(4, 0)] = 0;
        rem_step_by_32_M_128(rem, 28, term, q, rem);
        q = rem[index_word_hi(4)];

        if (q & 0x80000000) {
            sub_1_M_X96(extSigZ);
        } else {
            if (q || rem[index_word(4, 1)] || rem[index_word(4, 2)]) {
                extSigZ[index_word_lo(3)] |= 1;
            }
        }
    }

    round_pack_to_M_extF80M(0, expZ, extSigZ, extF80_rounding_precision, zPtr);
#endif
}
