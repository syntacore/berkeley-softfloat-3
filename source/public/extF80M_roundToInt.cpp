
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

#ifdef _MSC_VER
#pragma warning( disable : 5030)
#endif

void
extF80M_roundToInt(extFloat80_t const* const aPtr,
                   softfloat_round_mode roundingMode,
                   bool exact,
                   extFloat80_t* const zPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    *zPtr = extF80_roundToInt(*aPtr, roundingMode, exact);
#else
    using namespace softfloat::internals::slow_int64;

    uint16_t const signUI64 = static_cast<uint16_t>(aPtr->signExp & pack_to_extF80_UI64(1, 0));
    int32_t exp = get_exp(*aPtr);
    uint64_t sigA = aPtr->signif;

    if (0 == (sigA & INT64_MIN) && 0x7FFF != exp) {
        if (!sigA) {
            zPtr->signExp = signUI64;
            zPtr->signif = 0;
            return;
        }

        exp += norm_M_extF80Sig(&sigA);
    }

    if (exp <= 0x3FFE) {
        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        switch (roundingMode) {
        case softfloat_round_near_even:
            if (0 == (sigA & INT64_MAX)) {
                break;
            }

            [[fallthrough]];

        case softfloat_round_near_maxMag:
            if (0x3FFE == exp) {
                zPtr->signExp = uint16_t(signUI64 | 0x3FFF);
                zPtr->signif = static_cast<uint64_t>(INT64_MIN);
                return;
            }

            break;

        case softfloat_round_min:
            if (signUI64) {
                zPtr->signExp = static_cast<uint16_t>(signUI64 | UINT16_C(0x3FFF));
                zPtr->signif = static_cast<uint64_t>(INT64_MIN);
                return;
            }

            break;

        case softfloat_round_max:
            if (!signUI64) {
                zPtr->signExp = static_cast<uint16_t>(signUI64 | UINT16_C(0x3FFF));
                zPtr->signif = static_cast<uint64_t>(INT64_MIN);
                return;
            }

            break;
        }

        zPtr->signExp = signUI64;
        zPtr->signif = 0;
        return;
    }

    if (0x403E <= exp) {
        uint64_t sigZ;

        if (INT16_MAX == exp) {
            if (0 != (sigA & INT64_MAX)) {
                propagate_NaN_extF80M(aPtr, 0, zPtr);
                return;
            }

            sigZ = static_cast<uint64_t>(INT64_MIN);
        } else {
            sigZ = sigA;
        }

        /**
        @todo Warning   C4244   '=': conversion from 'int32_t' to 'uint16_t', possible loss of data
        */
        zPtr->signExp = static_cast<uint16_t>(signUI64 | exp);
        zPtr->signif = sigZ;
        return;
    }

    /**
    @todo Warning   C4244   '=': conversion from 'int32_t' to 'uint16_t', possible loss of data
    */
    uint16_t uiZ64 = static_cast<uint16_t>(signUI64 | exp);
    uint64_t const lastBitMask = static_cast<uint64_t>(1) << (0x403E - exp);
    uint64_t const roundBitsMask = lastBitMask - 1;
    uint64_t sigZ = sigA;

    if (softfloat_round_near_maxMag == roundingMode) {
        sigZ += lastBitMask >> 1;
    } else if (softfloat_round_near_even == roundingMode) {
        sigZ += lastBitMask >> 1;

        if (0 == (sigZ & roundBitsMask)) {
            sigZ &= ~lastBitMask;
        }
    } else if (softfloat_round_minMag != roundingMode && ((0 != signUI64) != (softfloat_round_max == roundingMode))) {
        sigZ += roundBitsMask;
    }

    sigZ &= ~roundBitsMask;

    if (!sigZ) {
        ++uiZ64;
        sigZ = UINT64_C(0x8000000000000000);
    }

    if (exact && (sigZ != sigA)) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    zPtr->signExp = uiZ64;
    zPtr->signif = sigZ;
    return;
#endif
}
