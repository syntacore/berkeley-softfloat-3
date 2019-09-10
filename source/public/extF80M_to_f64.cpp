
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

float64_t
extF80M_to_f64(extFloat80_t const* const aPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    return extF80_to_f64(*aPtr);
#else
    using namespace softfloat::internals::slow_int64;

    extFloat80M const* aSPtr = aPtr;
    bool const sign = is_sign(*aSPtr);
    int32_t exp = get_exp(*aSPtr);
    uint64_t sig = aSPtr->signif;

    if (0x7FFF == exp) {
        if (0 != (sig & UINT64_C(0x7FFFFFFFFFFFFFFF))) {
            return u_as_f(commonNaN_to_F64UI(commonNaN{*aSPtr}));
        }

        return u_as_f(pack_to_F64_UI(sign, 0x7FF, 0));
    }

    if (0 == (sig & UINT64_C(0x8000000000000000))) {
        if (0 == sig) {
            return u_as_f(pack_to_F64_UI(sign, 0, 0));
        }

        exp += norm_M_extF80Sig(&sig);
    }

    sig = short_shift_right_jam_64(sig, 1);
    exp -= 0x3C01;

    if (exp < -0x1000) {
        exp = -0x1000;
    }

    /**
    @todo Warning   C4242   'function': conversion from 'int32_t' to 'int16_t', possible loss of data
    */
    return round_pack_to_F64(sign, static_cast<int16_t>(exp), sig);
#endif
}
