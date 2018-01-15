
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

#include "primitives/functions.hpp"

uint32_t
softfloat_approxRecipSqrt32_1(uint32_t oddExpA,
                              uint32_t a)
{
    auto const index = (a >> 27 & 0xE) + oddExpA;
    uint16_t const eps = static_cast<uint16_t>(a >> 12);
    uint16_t const r0 =
        softfloat_approxRecipSqrt_1k0s[index] -
        ((softfloat_approxRecipSqrt_1k1s[index] * (uint32_t)eps) >> 20);
    uint32_t ESqrR0 = ((uint32_t)r0 * r0) << (!oddExpA ? 1 : 0);
    uint32_t const sigma0 = ~(uint32_t)(((uint32_t)ESqrR0 * (uint64_t)a) >> 23);
    uint32_t r =
        ((uint32_t)r0 << 16) +
        ((r0 * (uint64_t)sigma0) >> 25);
    uint32_t const sqrSigma0 = ((uint64_t)sigma0 * sigma0) >> 32;
    r += ((uint32_t)((r >> 1) + (r >> 3) - ((uint32_t)r0 << 14)) * (uint64_t)sqrSigma0) >> 48;
    return !(r & 0x80000000) ? 0x80000000 : r;
}
