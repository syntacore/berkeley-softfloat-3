
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

void
extF80M_mul(extFloat80_t const* const aPtr,
            extFloat80_t const* const bPtr,
            extFloat80_t* const zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    *zPtr = extF80_mul(*aPtr, *bPtr);
#else
    using namespace softfloat::internals::slow_int64;

    int32_t const expA = get_exp(*aPtr);
    int32_t const expB = get_exp(*bPtr);
    bool const signZ = is_sign(*aPtr) != is_sign(*bPtr);

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (!try_propagate_NaN_M_extF80(aPtr, bPtr, zPtr)) {
            if ((0 == aPtr->signif && 0x7FFF != expA) || (0 == bPtr->signif && 0x7FFF != expB)) {
                invalid_M_extF80(zPtr);
            } else {
                zPtr->signExp = pack_to_extF80_UI64(signZ, 0x7FFF);
                zPtr->signif = UINT64_C(0x8000000000000000);
            }
        }
    } else {
        auto expA_1 = 0 == expA ? 1 : expA;

        uint64_t sigA = aPtr->signif;

        if (0 == (sigA & UINT64_C(0x8000000000000000))) {
            if (0 == sigA) {
                zPtr->signExp = pack_to_extF80_UI64(signZ, 0);
                zPtr->signif = 0;
                return;
            }

            expA_1 += norm_M_extF80Sig(&sigA);
        }

        auto expB_1 = 0 == expB ? 1 : expB;

        uint64_t sigB = bPtr->signif;

        if (!(sigB & UINT64_C(0x8000000000000000))) {
            if (!sigB) {
                zPtr->signExp = pack_to_extF80_UI64(signZ, 0);
                zPtr->signif = 0;
                return;
            }

            expB_1 += norm_M_extF80Sig(&sigB);
        }

        int32_t const expZ = expA_1 + expB_1 - 0x3FFE;
        uint32_t sigProd[4];
        mul_M_64_to_128(sigA, sigB, sigProd);

        if (0 != sigProd[index_word_lo(4)]) {
            sigProd[index_word(4, 1)] |= 1;
        }

        auto const extSigZPtr = &sigProd[index_multiword_hi(4, 3)];

        if (sigProd[index_word_hi(4)] < 0x80000000) {
            add_M_96(extSigZPtr, extSigZPtr, extSigZPtr);
            round_pack_to_M_extF80M(signZ, expZ - 1, extSigZPtr, extF80_rounding_precision, zPtr);
        } else {
            round_pack_to_M_extF80M(signZ, expZ, extSigZPtr, extF80_rounding_precision, zPtr);
        }
    }
#endif
}
