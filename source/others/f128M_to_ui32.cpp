
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

uint32_t
f128M_to_ui32(float128_t const* const aPtr,
              softfloat_round_mode const roundingMode,
              bool const exact)
{
#ifdef SOFTFLOAT_FAST_INT64
    return f128_to_ui32(*aPtr, roundingMode, exact);
#else
    using namespace softfloat::internals::slow_int64;

    static bool const fromNaN_is_same_as_pos_overflow = ui32_fromNaN == ui32_fromPosOverflow;
    static bool const fromNaN_is_same_as_neg_overflow = ui32_fromNaN == ui32_fromNegOverflow;
    static bool const fromNaN_is_same_as_both_overflow = fromNaN_is_same_as_pos_overflow && fromNaN_is_same_as_neg_overflow;
    static bool const fromNaN_is_same_as_any_overflow = fromNaN_is_same_as_pos_overflow || fromNaN_is_same_as_neg_overflow;
    auto const aWPtr = reinterpret_cast<uint32_t const*>(aPtr);
    uint32_t const uiA96 = aWPtr[indexWordHi(4)];
    bool sign = is_sign(uiA96);
    int32_t const exp = expF128UI96(uiA96);
    uint64_t sig64 = static_cast<uint64_t>(fracF128UI96(uiA96)) << 32 | aWPtr[indexWord(4, 2)];

    if (aWPtr[indexWord(4, 1)] | aWPtr[indexWord(4, 0)]) {
        sig64 |= 1;
    }

    // TODO: check and re-factor
    if (!fromNaN_is_same_as_both_overflow && 0x7FFF == exp && 0 != sig64) {
        if (!fromNaN_is_same_as_any_overflow) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return ui32_fromNaN;
        }

        sign = fromNaN_is_same_as_neg_overflow;
    }

    if (exp) {
        sig64 |= UINT64_C(0x0001000000000000);
    }

    int32_t const shiftDist = 0x4023 - exp;
    return roundPackTo<uint32_t>(sign, 0 < shiftDist ? softfloat_shiftRightJam64(sig64, static_cast<uint32_t>(shiftDist)) : sig64, roundingMode, exact);
#endif
}
