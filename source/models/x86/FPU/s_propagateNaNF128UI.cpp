
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

namespace softfloat {
namespace internals {
namespace fast_int64 {
namespace x86 {
namespace {
static uint64_t const non_signal_bit = UINT64_C(0x0000'8000'0000'0000);
static uint64_t const suppress_sign_mask = UINT64_C(0x7FFF'FFFF'FFFF'FFFF);
}

uint128
propagate_NaN(uint64_t const& uiA64,
              uint64_t const& uiA0,
              uint64_t const& uiB64,
              uint64_t const& uiB0)
{
    bool const isSigNaNA = is_sNaN(uiA64, uiA0);
    bool const isSigNaNB = is_sNaN(uiB64, uiB0);
    /* Make NaNs non-signaling. */
    uint64_t const uiNonsigA64 = uiA64 | non_signal_bit;
    uint64_t const uiNonsigB64 = uiB64 | non_signal_bit;
    uint128 nsA(uiNonsigA64, uiA0);
    uint128 nsB(uiNonsigB64, uiB0);

    if (isSigNaNA || isSigNaNB) {
        softfloat_raiseFlags(softfloat_flag_invalid);

        if (!isSigNaNA) {
            return is_NaN(uiA64, uiA0) ? nsA : nsB;
        }

        assert(isSigNaNA);

        if (!isSigNaNB) {
            return is_NaN(uiB64, uiB0) ? nsB : nsA;
        }

        assert(isSigNaNA && isSigNaNB);
    }

    uint64_t const uiMagA64 = uiNonsigA64 & suppress_sign_mask;
    uint64_t const uiMagB64 = uiNonsigB64 & suppress_sign_mask;

    if (uiMagA64 < uiMagB64) {
        return nsB;
    }

    assert(uiMagB64 <= uiMagA64);

    if (uiMagB64 < uiMagA64) {
        return nsA;
    }

    assert(uiMagB64 == uiMagA64);

    if (uiA0 < uiB0) {
        return nsB;
    }

    assert(uiB0 <= uiA0);

    if (uiB0 < uiA0) {
        return nsA;
    }

    assert(uiB0 == uiA0);

    return
        uiNonsigA64 < uiNonsigB64 ? nsA : nsB;
}
}  // namespace x86
}  // namespace fast_int64
namespace slow_int64 {
namespace x86 {

namespace {
static inline void
result_copy(uint32_t* const zWPtr,
            uint32_t const* const ptr)
{
    zWPtr[indexWordHi(4)] = ptr[indexWordHi(4)] | 0x00008000;
    zWPtr[indexWord(4, 2)] = ptr[indexWord(4, 2)];
    zWPtr[indexWord(4, 1)] = ptr[indexWord(4, 1)];
    zWPtr[indexWord(4, 0)] = ptr[indexWord(4, 0)];
}
}  // namespace

/**
Assuming at least one of the two 128-bit floating-point values pointed to by
`aWPtr' and `bWPtr' is a NaN, stores the combined NaN result at the location
pointed to by `zWPtr'.  If either original floating-point value is a
signaling NaN, the invalid exception is raised.  Each of `aWPtr', `bWPtr',
and `zWPtr' points to an array of four 32-bit elements that concatenate in
the platform's normal endian order to form a 128-bit floating-point value.
*/
void
propagate_NaN_F128M(uint32_t const* const aWPtr,
                            uint32_t const* const bWPtr,
                            uint32_t* const zWPtr)
{

    bool const isSigNaNA = f128M_isSignalingNaN(reinterpret_cast<float128_t const*>(aWPtr));
    bool const isSigNaNB = f128M_isSignalingNaN(reinterpret_cast<float128_t const*>(bWPtr));

    if (!bWPtr) {
        if (isSigNaNA) {
            softfloat_raiseFlags(softfloat_flag_invalid);
        }

        result_copy(zWPtr, aWPtr);
        return;
    }

    if (isSigNaNA || isSigNaNB) {
        softfloat_raiseFlags(softfloat_flag_invalid);

        if (!isSigNaNA) {
            result_copy(zWPtr, softfloat_isNaNF128M(aWPtr) ? aWPtr : bWPtr);
            return;
        }

        if (!isSigNaNB) {
            result_copy(zWPtr, softfloat_isNaNF128M(bWPtr) ? bWPtr : aWPtr);
            return;
        }
    }

    auto const uiA96 = aWPtr[indexWordHi(4)];
    auto const uiB96 = bWPtr[indexWordHi(4)];
    auto const wordMagA_3 = uiA96 & 0x7FFFFFFF;
    auto const wordMagB_3 = uiB96 & 0x7FFFFFFF;
    auto const wordMagA_2 = aWPtr[indexWord(4, 2)];
    auto const wordMagB_2 = bWPtr[indexWord(4, 2)];
    auto const wordMagA_1 = aWPtr[indexWord(4, 1)];
    auto const wordMagB_1 = bWPtr[indexWord(4, 1)];
    auto const wordMagA_0 = aWPtr[indexWord(4, 0)];
    auto const wordMagB_0 = bWPtr[indexWord(4, 0)];
    bool const isA =
        (!(wordMagA_3 < wordMagB_3) && wordMagB_3 < wordMagA_3) ||
        (
            !(wordMagA_2 < wordMagB_2) &&
            (
                wordMagB_2 < wordMagA_2 ||
                (
                    !(wordMagA_1 < wordMagB_1) &&
                    (
                        wordMagB_1 < wordMagA_1 ||
                        (
                            !(wordMagA_0 < wordMagB_0) &&
                            (wordMagB_0 < wordMagA_0 || uiA96 < uiB96)
                        )
                    )
                )

            )
        );
    result_copy(zWPtr, isA ? aWPtr : bWPtr);
}

}  // namespace x86
}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
