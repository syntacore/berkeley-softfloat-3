
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

int32_t
extF80M_to_i32(extFloat80_t const* const aPtr,
               softfloat_round_mode const roundingMode,
               bool const exact)
{
#if (SOFTFLOAT_FAST_INT64)
    return extF80_to_i32(*aPtr, roundingMode, exact);
#else
    using namespace softfloat::internals;

    bool const sign = is_sign(*aPtr);
    int32_t const exp = get_exp(*aPtr);
    uint64_t sig = aPtr->signif;
    int32_t const shiftDist = 0x4032 - exp;

    if (shiftDist <= 0) {
        if (sig >> 32) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return
                0x7FFF == exp && 0 != (sig & UINT64_C(0x7FFFFFFFFFFFFFFF)) ? i32_fromNaN :
                sign ? i32_fromNegOverflow :
                i32_fromPosOverflow;
        }

        if (-32 < shiftDist) {
            sig <<= -shiftDist;
        } else if (0 != sig) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return
                0x7FFF == exp && 0 != (sig & UINT64_C(0x7FFFFFFFFFFFFFFF)) ? i32_fromNaN :
                sign ? i32_fromNegOverflow :
                i32_fromPosOverflow;
        }
    } else {
        sig = shift_right_jam_64(sig, static_cast<uint32_t>(shiftDist));
    }

    return round_pack_to<int32_t>(sign, sig, roundingMode, exact);
#endif
}
