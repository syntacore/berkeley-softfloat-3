
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

uint32_t
extF80M_to_ui32_r_minMag(extFloat80_t const* const aPtr,
                         bool const exact)
{
    using namespace softfloat::internals;
    uint16_t uiA64 = aPtr->signExp;
    int32_t exp = expExtF80UI64(uiA64);
    uint64_t sig = aPtr->signif;

    if (!sig && exp != INT16_MAX) {
        return 0;
    }

    int32_t const shiftDist = 0x403E - exp;

    if (64 <= shiftDist) {
        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return 0;
    }

    bool const sign = is_sign(uiA64);

    if (shiftDist < 0) {
        if (sign || sig >> 32 || (shiftDist <= -31)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return
                exp == INT16_MAX && 0 != (sig & INT64_MAX) ? ui32_fromNaN :
                sign ? ui32_fromNegOverflow : ui32_fromPosOverflow;
        }

        uint64_t const shiftedSig = static_cast<uint64_t>(static_cast<uint32_t>(sig)) << -shiftDist;

        if (0 != (shiftedSig >> 32)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return
                exp == INT16_MAX && 0 != (sig & INT64_MAX) ? ui32_fromNaN :
                sign ? ui32_fromNegOverflow : ui32_fromPosOverflow;
        }

        return static_cast<uint32_t>(shiftedSig);
    }

    uint64_t const shiftedSig = 0 != shiftDist ? sig >>= shiftDist : sig;

    if (0 != (shiftedSig >> 32)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            exp == INT16_MAX && 0 != (sig & INT64_MAX) ? ui32_fromNaN :
            sign ? ui32_fromNegOverflow : ui32_fromPosOverflow;
    }

    uint32_t const z = static_cast<uint32_t>(shiftedSig);

    if (sign && 0 != z) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return
            exp == INT16_MAX && 0 != (sig & INT64_MAX) ? ui32_fromNaN :
            sign ? ui32_fromNegOverflow : ui32_fromPosOverflow;
    }

    if (exact && shiftDist && (static_cast<uint64_t>(z) << shiftDist != sig)) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return z;
}
