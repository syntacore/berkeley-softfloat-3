
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
extF80M_to_ui64_r_minMag(extFloat80_t const* const aPtr,
                         bool const exact)
{
#if (SOFTFLOAT_FAST_INT64)
    return extF80_to_ui64_r_minMag(*aPtr, exact);
#else
    using namespace softfloat::internals;

    extFloat80M const* const aSPtr = aPtr;
    int32_t const exp = get_exp(*aPtr);

    if (0 == aSPtr->signif && 0x7FFF != exp) {
        return 0;
    }

    int32_t const shiftDist = 0x403E - exp;

    if (64 <= shiftDist) {
        if (exact) {
            softfloat_raiseFlags(softfloat_flag_inexact);
        }

        return 0;
    }

    bool const sign = is_sign(*aSPtr);

    if (shiftDist < 0) {
        if (!sign && -63 < shiftDist) {
            auto const shiftDist_1 = -shiftDist;
            uint64_t const z = aSPtr->signif << shiftDist_1;

            if (z >> shiftDist_1 == aSPtr->signif) {
                return z;
            }
        }
    } else {
        assert(0 <= shiftDist && shiftDist < 64);
        uint64_t const z = aSPtr->signif >> shiftDist;

        if (!sign || 0 == z) {
            if (exact && 0 != shiftDist && z << shiftDist != aSPtr->signif) {
                softfloat_raiseFlags(softfloat_flag_inexact);
            }

            return z;
        }
    }

    softfloat_raiseFlags(softfloat_flag_invalid);
    return
        0x7FFF == exp && 0 != (aSPtr->signif & UINT64_C(0x7FFFFFFFFFFFFFFF)) ?
        ui64_fromNaN :
        sign ?
        ui64_fromNegOverflow :
        ui64_fromPosOverflow;
#endif
}
