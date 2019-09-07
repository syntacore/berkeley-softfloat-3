
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

uint64_t
f64_to_ui64(float64_t const a,
            softfloat_round_mode const roundingMode,
            bool const exact)
{
#ifdef SOFTFLOAT_FAST_INT64
    using namespace softfloat::internals::fast_int64;

    bool const sign = is_sign(a);
    int16_t const exp = get_exp(a);
    uint64_t const sig = get_frac(a) | (0 != exp ? UINT64_C(0x0010000000000000) : UINT64_C(0));
    int16_t const shiftDist = 0x433 - exp;

    if (shiftDist <= 0) {
        if (shiftDist < -11) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return
                is_NaN(a) ? ui64_fromNaN :
                sign ? ui64_fromNegOverflow : ui64_fromPosOverflow;
        }

        return roundPackTo<uint64_t>(sign, sig << -shiftDist, 0, roundingMode, exact);
    }

    uint64_extra const sigExtra = softfloat_shiftRightJam64Extra(sig, 0, static_cast<uint32_t>(shiftDist));
    return roundPackTo<uint64_t>(sign, sigExtra.v, sigExtra.extra, roundingMode, exact);
#else
    using namespace softfloat::internals::slow_int64;

    bool const sign = is_sign(a);
    int16_t const exp = get_exp(a);
    uint64_t const sig = get_frac(a) | (0 != exp ? UINT64_C(0x0010000000000000) : UINT64_C(0));
    int16_t const shiftDist = 0x433 - exp;

    if (shiftDist <= 0) {
        if (shiftDist < -11) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return
                is_NaN(a) ? ui64_fromNaN :
                sign ? ui64_fromNegOverflow : ui64_fromPosOverflow;
        }

        auto const sig_1 = sig << -shiftDist;
        uint32_t extSig[3];
        extSig[indexWord(3, 0)] = 0;
        extSig[indexWord(3, 2)] = sig_1 >> 32;
        extSig[indexWord(3, 1)] = static_cast<uint32_t>(sig_1);
        return roundPackMTo<uint64_t>(sign, extSig, roundingMode, exact);
    }

    uint32_t extSig[3];
    extSig[indexWord(3, 0)] = 0;
    extSig[indexWord(3, 2)] = sig >> 32;
    extSig[indexWord(3, 1)] = static_cast<uint32_t>(sig);
    softfloat_shiftRightJam96M(extSig, static_cast<uint8_t>(shiftDist), extSig);
    return roundPackMTo<uint64_t>(sign, extSig, roundingMode, exact);
#endif
}
