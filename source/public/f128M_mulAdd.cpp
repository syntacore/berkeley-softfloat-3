
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014 The Regents of the University of California.
All rights reserved.

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

void
f128M_mulAdd(float128_t const* const aPtr,
             float128_t const* const bPtr,
             float128_t const* const cPtr,
             float128_t* const zPtr)
{
#ifdef SOFTFLOAT_FAST_INT64
    using namespace softfloat::internals::fast_int64;

    uint64_t const* const aWPtr = reinterpret_cast<uint64_t const*>(aPtr);
    uint64_t const* const bWPtr = reinterpret_cast<uint64_t const*>(bPtr);
    uint64_t const* const cWPtr = reinterpret_cast<uint64_t const*>(cPtr);
    uint64_t const uiA64 = aWPtr[indexWord(2, 1)];
    uint64_t const uiA0 = aWPtr[indexWord(2, 0)];
    uint64_t const uiB64 = bWPtr[indexWord(2, 1)];
    uint64_t const uiB0 = bWPtr[indexWord(2, 0)];
    uint64_t const uiC64 = cWPtr[indexWord(2, 1)];
    uint64_t const uiC0 = cWPtr[indexWord(2, 0)];
    *zPtr = softfloat_mulAddF128(softfloat_mulAdd_madd,
                                 uiA64, uiA0,
                                 uiB64, uiB0,
                                 uiC64, uiC0);

#else
    using namespace softfloat::internals::slow_int64;

    softfloat_mulAddF128M(softfloat_mulAdd_madd,
                          reinterpret_cast<const uint32_t*>(aPtr),
                          reinterpret_cast<const uint32_t*>(bPtr), 
                          reinterpret_cast<const uint32_t*>(cPtr),
                          reinterpret_cast<uint32_t*>(zPtr));
#endif

}
