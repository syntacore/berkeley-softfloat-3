
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

#if (SOFTFLOAT_FAST_INT64)

void
extF80M_add(extFloat80_t const *const aPtr,
            extFloat80_t const *const bPtr,
            extFloat80_t *const zPtr)
{
    using namespace softfloat::internals::fast_int64;

    auto const uiA64 = aPtr->signExp;
    auto const uiA0 = aPtr->signif;
    bool const signA = is_sign(uiA64);
    auto const uiB64 = bPtr->signExp;
    auto const uiB0 = bPtr->signif;
    bool const signB = is_sign(uiB64);
    *zPtr =
        signA == signB ? add_magnitudes(uiA64, uiA0, uiB64, uiB0, signA) :
        sub_magnitudes(uiA64, uiA0, uiB64, uiB0, signA);
}

#else

void
extF80M_add(extFloat80_t const *const aPtr,
            extFloat80_t const *const bPtr,
            extFloat80_t *const zPtr)
{
    using namespace softfloat::internals::slow_int64;
    add_M_extF80(aPtr, bPtr, zPtr, false);
}

#endif
