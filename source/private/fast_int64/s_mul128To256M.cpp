
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015 The Regents of the University of
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

namespace softfloat {
namespace internals {
namespace fast_int64 {

void
mul_M_128_to_256(uint64_t const& a64,
             uint64_t const& a0,
             uint64_t const& b64,
             uint64_t const& b0,
             uint64_t* zPtr)
{
    uint128 const p0 = mul_64_to_128(a0, b0);
    zPtr[index_word(4, 0)] = p0.v0;
    uint128 p64 = mul_64_to_128(a64, b0);
    uint64_t z64 = p64.v0 + p0.v64;
    uint64_t z128 = p64.v64 + (z64 < p64.v0);
    uint128 const p128 = mul_64_to_128(a64, b64);
    z128 += p128.v0;
    uint64_t const z192 = p128.v64 + (z128 < p128.v0);
    p64 = mul_64_to_128(a0, b64);
    z64 += p64.v0;
    zPtr[index_word(4, 1)] = z64;
    p64.v64 += (z64 < p64.v0);
    z128 += p64.v64;
    zPtr[index_word(4, 2)] = z128;
    zPtr[index_word(4, 3)] = z192 + (z128 < p64.v64);
}

}  // namespace fast_int64
}  // namespace internals
}  // namespace softfloat
