
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

int64_t
f128M_to_i64(float128_t const* aPtr,
             softfloat_round_mode const roundingMode,
             bool const exact)
{
#if (SOFTFLOAT_FAST_INT64)
    return f128_to_i64(*aPtr, roundingMode, exact);
#else
    using namespace softfloat::internals::slow_int64;

    auto const aWPtr = reinterpret_cast<uint32_t const*>(aPtr);
    uint32_t const uiA96 = aWPtr[index_word_hi(4)];
    bool const sign = is_sign(uiA96);
    int32_t const exp = exp_F128_UI96(uiA96);
    uint32_t sig96 = frac_F128_UI96(uiA96);
    int32_t const shiftDist = 0x404F - exp;

    if (shiftDist < 17) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            INT32_C(0x7FFF) == exp && 0 != (sig96 || (aWPtr[index_word(4, 2)] | aWPtr[index_word(4, 1)] | aWPtr[index_word(4, 0)])) ? i64_fromNaN :
            sign ? i64_fromNegOverflow :
            i64_fromPosOverflow;
    }

    if (exp) {
        sig96 |= UINT32_C(0x00010000);
    }

    uint32_t sig[4];
    sig[index_word(4, 3)] = sig96;
    sig[index_word(4, 2)] = aWPtr[index_word(4, 2)];
    sig[index_word(4, 1)] = aWPtr[index_word(4, 1)];
    sig[index_word(4, 0)] = aWPtr[index_word(4, 0)];
    shift_right_jam_M_128(sig, static_cast<uint8_t>(shiftDist), sig);
    return
        round_pack_to_M<int64_t>(sign, sig + index_multiword_lo(4, 3), roundingMode, exact);

#endif
}
