
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

#ifdef SOFTFLOAT_FAST_INT64
int32_t
extF80M_to_i32_r_minMag(extFloat80_t const* const aPtr,
                        bool exact)
{
    return extF80_to_i32_r_minMag(*aPtr, exact);
}

#else

namespace {
static int32_t
invalid(bool sign, int32_t exp, uint64_t sig)
{
    using namespace softfloat::internals;
    softfloat_raiseFlags(softfloat_flag_invalid);
    return
        exp == INT16_MAX && (sig & INT64_MAX) ?
        i32_fromNaN :
        sign ?
        i32_fromNegOverflow :
        i32_fromPosOverflow;
}
}  // namespace

int32_t
extF80M_to_i32_r_minMag(extFloat80_t const* const aPtr,
                        bool exact)
{
    using namespace softfloat::internals;
    uint16_t const uiA64 = aPtr->signExp;
    int32_t const exp = expExtF80UI64(uiA64);
    uint64_t const sig = aPtr->signif;

    if (0 == sig && 0x7FFF != exp) {
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
    bool raiseInexact = false;

    uint32_t absZ;

    if (shiftDist < 0) {
        if (sig >> 32 || (shiftDist <= -31)) {
            return invalid(sign, exp, sig);
        }

        uint64_t const shiftedSig = static_cast<uint64_t>(static_cast<uint32_t>(sig)) << -shiftDist;

        if (shiftedSig >> 32) {
            return invalid(sign, exp, sig);
        }

        absZ = static_cast<uint32_t>(shiftedSig);
    } else {
        uint64_t const shiftedSig = sig >> shiftDist;

        if (0 != (shiftedSig >> 32)) {
            return invalid(sign, exp, sig);
        }

        absZ = static_cast<uint32_t>(shiftedSig);

        if (exact && 0 != shiftDist) {
            raiseInexact = static_cast<uint64_t>(absZ) << shiftDist != sig;
        }
    }

    if (sign) {
        if (0x80000000 < absZ) {
            return invalid(sign, exp, sig);
        }

        if (raiseInexact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return -static_cast<int32_t>(absZ);
    }

    if (0x80000000u <= absZ) {
        return invalid(sign, exp, sig);
    }

    if (raiseInexact) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return static_cast<int32_t>(absZ);
}
#endif
