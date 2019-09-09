
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
All rights reserved.
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

bool
extF80M_le_quiet(extFloat80_t const* const aPtr,
                 extFloat80_t const* const bPtr)
{
#if (SOFTFLOAT_FAST_INT64)
    return extF80_le_quiet(*aPtr, *bPtr);
#else
    using namespace softfloat::internals::slow_int64;

    if (is_NaN(*aPtr) || is_NaN(*bPtr)) {
        if (is_sNaN(*aPtr) || is_sNaN(*bPtr)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
        }

        return false;
    }

    uint16_t const uiA64 = aPtr->signExp;
    uint64_t const uiA0 = aPtr->signif;
    uint16_t const uiB64 = bPtr->signExp;
    uint64_t const uiB0 = bPtr->signif;

    bool const signA = is_sign(uiA64);

    if (0 != ((uiA64 ^ uiB64) & 0x8000)) {
        /* Signs are different. */
        return signA || 0 == (uiA0 | uiB0);
    }

    if (0 == ((uiA0 & uiB0) & UINT64_C(0x8000000000000000))) {
        return (compare_non_norm_M_extF80(aPtr, bPtr) <= 0);
    }

    if (uiA64 == uiB64) {
        return uiA0 == uiB0 || (signA != (uiA0 < uiB0));
    }

    return signA != (uiA64 < uiB64);
#endif
}
