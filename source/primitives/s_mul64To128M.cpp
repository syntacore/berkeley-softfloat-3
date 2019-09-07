
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
namespace slow_int64 {

void
softfloat_mul64To128M(uint64_t a,
                      uint64_t b,
                      uint32_t *zPtr)
{
    uint32_t const a32 = a >> 32;
    uint32_t const a0 = static_cast<uint32_t>(a);
    uint32_t const b32 = b >> 32;
    uint32_t const b0 = static_cast<uint32_t>(b);
    uint64_t const mid1 = static_cast<uint64_t>(a32) * b0;
    uint64_t const mid = mid1 + static_cast<uint64_t>(a0) * b32;
    uint64_t const mid2 = mid << 32;
    uint64_t const z0 = static_cast<uint64_t>(a0) * b0 + mid2;
    uint64_t const z64 = (static_cast<uint64_t>(a32) * b32 + ((static_cast<uint64_t>(mid < mid1) << 32) | (mid >> 32))) + (z0 < mid2);
    zPtr[indexWord(4, 1)] = z0 >> 32;
    zPtr[indexWord(4, 0)] = static_cast<uint32_t>(z0);
    zPtr[indexWord(4, 3)] = z64 >> 32;
    zPtr[indexWord(4, 2)] = static_cast<uint32_t>(z64);
}

}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
