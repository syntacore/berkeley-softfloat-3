
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

#if !(SOFTFLOAT_FAST_INT64)
namespace {
static uint32_t
invalid(bool const sign,
        int32_t const exp,
        uint64_t const sig)
{
    using namespace softfloat::internals;

    return
        INT16_MAX == exp && 0 != (sig & INT64_MAX) ?
        ui32_fromNaN :
        sign ?
        ui32_fromNegOverflow :
        ui32_fromPosOverflow;
}
}  // namespace
#endif

uint32_t
extF80M_to_ui32_r_minMag(extFloat80_t const* const aPtr,
                         bool const exact)
{
#if (SOFTFLOAT_FAST_INT64)
    return extF80_to_ui32_r_minMag(*aPtr, exact);
#else
    using namespace softfloat::internals;

    int32_t const exp = exp_extF80_UI64(aPtr->signExp);
    uint64_t const sig = aPtr->signif;

    if (0 == sig && INT16_MAX != exp) {
        return 0;
    }

    int32_t const shiftDist = 0x403E - exp;

    if (64 <= shiftDist) {
        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return 0;
    }

    bool const sign = is_sign(aPtr->signExp);

    if (shiftDist < 0) {
        if (sign || sig >> 32 || (shiftDist <= -31)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return invalid(sign, exp, sig);
        }

        uint64_t const shiftedSig = static_cast<uint64_t>(static_cast<uint32_t>(sig)) << -shiftDist;

        if (0 != (shiftedSig >> 32)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return invalid(sign, exp, sig);
        }

        return static_cast<uint32_t>(shiftedSig);
    }

    uint64_t const shiftedSig = 0 != shiftDist ? sig >> shiftDist : sig;

    if (0 != (shiftedSig >> 32)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return invalid(sign, exp, sig);
    }

    uint32_t const z = static_cast<uint32_t>(shiftedSig);

    if (0 != sign && 0 != z) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return invalid(sign, exp, sig);
    }

    if (exact && 0 != shiftDist && (static_cast<uint64_t>(z) << shiftDist != sig)) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return z;
#endif
}
