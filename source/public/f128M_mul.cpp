
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

void
f128M_mul(float128_t const* const aPtr,
          float128_t const* const bPtr,
          float128_t* const zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    *zPtr = f128_mul(*aPtr, *bPtr);
#else
    using namespace softfloat::internals::slow_int64;

    uint32_t const* const aWPtr = reinterpret_cast<const uint32_t*>(aPtr);
    uint32_t const* const bWPtr = reinterpret_cast<const uint32_t*>(bPtr);
    uint32_t* const zWPtr = reinterpret_cast<uint32_t*>(zPtr);

    uint32_t const uiA96 = aWPtr[index_word_hi(4)];
    int32_t expA = exp_F128_UI96(uiA96);
    uint32_t const uiB96 = bWPtr[index_word_hi(4)];
    int32_t expB = exp_F128_UI96(uiB96);
    bool const signZ = is_sign(uiA96) != is_sign(uiB96);

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (try_propagate_NaN_M_F128(aWPtr, bWPtr, zWPtr)) {
            return;
        }

        if (0 == expA || 0 == expB) {
            uint32_t const* ptr = 0 == expA ? aWPtr : bWPtr;

            if (
                0 == frac_F128_UI96(ptr[index_word_hi(4)]) &&
                0 == (ptr[index_word(4, 2)] | ptr[index_word(4, 1)] | ptr[index_word(4, 0)])
            ) {
                invalid_M_F128(zWPtr);
                return;
            }
        }

        zWPtr[index_word_hi(4)] = pack_to_F128_UI96(signZ, 0x7FFF, 0);
        zWPtr[index_word(4, 2)] = 0;
        zWPtr[index_word(4, 1)] = 0;
        zWPtr[index_word(4, 0)] = 0;
        return;
    }

    uint32_t sigA[4];

    if (0 != expA) {
        sigA[index_word_hi(4)] = frac_F128_UI96(uiA96) | 0x00010000;
        sigA[index_word(4, 2)] = aWPtr[index_word(4, 2)];
        sigA[index_word(4, 1)] = aWPtr[index_word(4, 1)];
        sigA[index_word(4, 0)] = aWPtr[index_word(4, 0)];
    } else {
        expA = shift_norm_sig_M_F128(aWPtr, 0, sigA);

        if (expA == -128) {
            zWPtr[index_word_hi(4)] = pack_to_F128_UI96(signZ, 0, 0);
            zWPtr[index_word(4, 2)] = 0;
            zWPtr[index_word(4, 1)] = 0;
            zWPtr[index_word(4, 0)] = 0;
            return;
        }
    }

    uint32_t sigB[4];

    if (0 != expB) {
        sigB[index_word_hi(4)] = frac_F128_UI96(uiB96) | 0x00010000;
        sigB[index_word(4, 2)] = bWPtr[index_word(4, 2)];
        sigB[index_word(4, 1)] = bWPtr[index_word(4, 1)];
        sigB[index_word(4, 0)] = bWPtr[index_word(4, 0)];
    } else {
        expB = shift_norm_sig_M_F128(bWPtr, 0, sigB);

        if (expB == -128) {
            zWPtr[index_word_hi(4)] = pack_to_F128_UI96(signZ, 0, 0);
            zWPtr[index_word(4, 2)] = 0;
            zWPtr[index_word(4, 1)] = 0;
            zWPtr[index_word(4, 0)] = 0;
            return;
        }
    }

    int32_t expZ = expA + expB - 0x4000;
    uint32_t sigProd[8];
    mul_M_128_to_256(sigA, sigB, sigProd);

    if (0 != sigProd[index_word(8, 2)] || 0 != (sigProd[index_word(8, 1)] | sigProd[index_word(8, 0)])) {
        sigProd[index_word(8, 3)] |= 1;
    }

    uint32_t* const extSigZPtr = &sigProd[index_multiword_hi(8, 5)];
    uint8_t shiftDist = 16;

    if (0 != (extSigZPtr[index_word_hi(5)] & 2)) {
        ++expZ;
        shiftDist = 15;
    }

    short_shift_left_M_160(extSigZPtr, shiftDist, extSigZPtr);
    round_pack_to_M_F128(signZ, expZ, extSigZPtr, zWPtr);
#endif
}
