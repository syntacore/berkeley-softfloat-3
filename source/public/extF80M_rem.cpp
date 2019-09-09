
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
extF80M_rem(extFloat80_t const* const aPtr,
            extFloat80_t const* const bPtr,
            extFloat80_t* const zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    *zPtr = extF80_rem(*aPtr, *bPtr);
#else
    using namespace softfloat::internals::slow_int64;
    uint16_t const uiA64 = aPtr->signExp;
    int32_t expA = exp_extF80_UI64(uiA64);
    int32_t expB = exp_extF80_UI64(bPtr->signExp);

    if (0x7FFF == expA || 0x7FFF == expB) {
        if (try_propagate_NaN_M_extF80(aPtr, bPtr, zPtr)) {
            return;
        }

        if (0x7FFF == expA) {
            invalid_M_extF80(zPtr);
            return;
        }

        /*
        If we get here, then argument b is an infinity and `expB' is 0x7FFF;
        Doubling `expB' is an easy way to ensure that `expDiff' later is
        less than -1, which will result in returning a canonicalized version
        of argument a.
        */
        expB += expB;
    }

    if (0 == expB) {
        expB = 1;
    }

    uint64_t x64 = bPtr->signif;

    if (0 == (x64 & UINT64_C(0x8000000000000000))) {
        if (0 == x64) {
            invalid_M_extF80(zPtr);
            return;
        }

        expB += norm_M_extF80Sig(&x64);
    }

    bool signRem = is_sign(uiA64);

    if (0 == expA) {
        expA = 1;
    }

    uint64_t sigA = aPtr->signif;

    if (0 == (sigA & UINT64_C(0x8000000000000000))) {
        if (0 == sigA) {
            zPtr->signExp = pack_to_extF80_UI64(signRem, 0);
            zPtr->signif = sigA >> 1;
            return;
        }

        expA += norm_M_extF80Sig(&sigA);
    }

    int32_t expDiff = expA - expB;

    if (expDiff < -1) {
        if (expA < 1) {
            sigA >>= 1 - expA;
            expA = 0;
        }

        zPtr->signExp = pack_to_extF80_UI64(signRem, static_cast<uint16_t>(expA));
        zPtr->signif = sigA;
    } else {
        uint32_t rem[3];
        rem[index_word(3, 2)] = sigA >> 34;
        rem[index_word(3, 1)] = static_cast<uint32_t>(sigA >> 2);
        rem[index_word(3, 0)] = static_cast<uint32_t>(sigA) << 30;

        uint32_t x[3];
        x[index_word(3, 0)] = static_cast<uint32_t>(x64) << 30;
        uint32_t const sig32B = x64 >> 32;
        x64 >>= 2;
        x[index_word(3, 2)] = x64 >> 32;
        x[index_word(3, 1)] = static_cast<uint32_t>(x64);

        uint32_t q;

        if (expDiff < 1) {
            if (expDiff) {
                --expB;
                add_M_96(x, x, x);
                q = 0;
            } else {
                q = 0u + !!(compare_M_96(x, rem) <= 0);

                if (q) {
                    sub_M_96(rem, x, rem);
                }
            }
        } else {
            uint32_t const recip32 = approx_recip_32_1(sig32B);
            expDiff -= 30;

            for (;;) {
                x64 = static_cast<uint64_t>(rem[index_word_hi(3)]) * recip32;

                if (expDiff < 0) {
                    break;
                } else {
                    q = (x64 + 0x80000000) >> 32;
                    rem_step_by_32_M_96(rem, 29, x, q, rem);

                    if (rem[index_word_hi(3)] & 0x80000000) {
                        add_M_96(rem, x, rem);
                    }

                    expDiff -= 29;
                }
            }

            /* `expDiff' cannot be less than -29 here. */
            q = static_cast<uint32_t>(x64 >> 32) >> (~expDiff & 31);
            rem_step_by_32_M_96(rem, static_cast<uint8_t>(expDiff + 30), x, q, rem);

            if (rem[index_word_hi(3)] & 0x80000000) {
                uint32_t* remPtr = rem;
                uint32_t rem2[3];
                uint32_t* altRemPtr = rem2;
                add_M_96(remPtr, x, altRemPtr);
                add_M_96(remPtr, altRemPtr, x);
                uint32_t const wordMeanRem = x[index_word_hi(3)];

                if (0 != (wordMeanRem & 0x80000000) || (0 == wordMeanRem && 0 != (q & 1) && !x[index_word(3, 0)] && !x[index_word(3, 1)])) {
                    remPtr = altRemPtr;
                }

                if (remPtr[index_word_hi(3)] & 0x80000000) {
                    signRem = !signRem;
                    neg_M_X96(remPtr);
                }

                norm_round_pack_to_M_extF80(signRem, expB + 2, remPtr, 80, zPtr);
                return;
            }
        }

        uint32_t* remPtr = rem;
        uint32_t rem2[3];
        uint32_t* altRemPtr = rem2;

        do {
            ++q;
            uint32_t* newRemPtr = altRemPtr;
            sub_M_96(remPtr, x, newRemPtr);
            altRemPtr = remPtr;
            remPtr = newRemPtr;
        } while (!(remPtr[index_word_hi(3)] & 0x80000000));

        add_M_96(remPtr, altRemPtr, x);
        uint32_t const wordMeanRem = x[index_word_hi(3)];

        if (0 != (wordMeanRem & 0x80000000) || (0 == wordMeanRem && 0 != (q & 1) && !x[index_word(3, 0)] && !x[index_word(3, 1)])) {
            remPtr = altRemPtr;
        }

        if (remPtr[index_word_hi(3)] & 0x80000000) {
            signRem = !signRem;
            neg_M_X96(remPtr);
        }

        norm_round_pack_to_M_extF80(signRem, expB + 2, remPtr, 80, zPtr);
    }
#endif
}
