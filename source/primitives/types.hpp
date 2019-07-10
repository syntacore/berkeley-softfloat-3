
/** @file

This C header file is part of the SoftFloat IEEE Floating-Point Arithmetic
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

#ifndef SOFTFLOAT_PRIMITIVES_TYPES_HPP_
#define SOFTFLOAT_PRIMITIVES_TYPES_HPP_

#include "softfloat/types.h"
#include <cstdint>
#include <climits>

namespace softfloat {
namespace internals {

/**
These macros are used to isolate the differences in word order between big-
endian and little-endian platforms.
*/
#ifdef BIG_ENDIAN

static int const wordIncr = -1;

static inline constexpr size_t
indexWord(size_t const total, size_t const n)
{
    return total - 1 - n;
}

static inline constexpr size_t
indexWordHi(size_t const total)
{
    return 0;
}

static inline constexpr size_t
indexWordLo(size_t const total)
{
    return total - 1;
}

static inline constexpr size_t
indexMultiword(size_t total,
               size_t m,
               size_t)
{
    return total - 1 - m;
}

static inline constexpr size_t
indexMultiwordHi(size_t,
                 size_t)
{
    return 0;
}

static inline constexpr size_t
indexMultiwordLo(size_t total,
                 size_t n)
{ 
    return total - n;
}

static inline constexpr size_t
indexMultiwordHiBut(size_t,
                    size_t) 
{
    return 0;
}

static inline constexpr size_t
indexMultiwordLoBut(size_t,
                    size_t n)
{
    return n;
}

#define INIT_UINTM4( v3, v2, v1, v0 ) { v3, v2, v1, v0 }

#else

static int const wordIncr = 1;

static inline constexpr size_t
indexWord(size_t, size_t const n)
{
    return n;
}

static inline constexpr size_t
indexWordHi(size_t const total)
{
    return total - 1u;
}

static inline constexpr size_t
indexWordLo(size_t const)
{
    return 0;
}

static inline constexpr size_t
indexMultiword(size_t,
               size_t,
               size_t n)
{
    return n;
}

static inline constexpr size_t
indexMultiwordHi(size_t total,
                 size_t n) 
{
    return total - n;
}

static inline constexpr size_t
indexMultiwordLo(size_t, size_t)
{
    return 0u;
}

static inline constexpr size_t
indexMultiwordHiBut(size_t,
                    size_t n)
{
    return n;
}

static inline constexpr size_t
indexMultiwordLoBut(size_t,
                    size_t)
{
    return 0u;
}

#define INIT_UINTM4( v3, v2, v1, v0 ) { v0, v1, v2, v3 }

#endif

static_assert(16 == CHAR_BIT * sizeof(float16_t), "Bad size of float16_t");
static_assert(32 == CHAR_BIT * sizeof(float32_t), "Bad size of float32_t");
static_assert(64 == CHAR_BIT * sizeof(float64_t), "Bad size of float64_t");
static_assert(80 <= CHAR_BIT * sizeof(extFloat80M) && CHAR_BIT * sizeof(extFloat80M) <= 128, "Bad size of extFloat80M");
static_assert(128 == CHAR_BIT * sizeof(float128_t), "Bad size of float128_t");

#ifdef SOFTFLOAT_FAST_INT64

struct uint128
{
    uint128() = default;

    constexpr
    uint128(uint64_t const& a64,
            uint64_t const& a0)
        : v64(a64)
        , v0(a0)
    {}

    explicit constexpr
    uint128(float128_t const& a)
        : v64(a.v64)
        , v0(a.v0)
    {}

    explicit constexpr
    operator float128_t()const
    {
#ifdef BIG_ENDIAN
        return float128_t{v64, v0};
#else
        return float128_t{v0, v64};
#endif
    }

    uint64_t v64;
    uint64_t v0;
};

static_assert(sizeof(uint128) == sizeof(float128_t), "sizeof(uint128) != sizeof(float128_t)");

struct uint64_extra
{
    uint64_extra() = default;

    constexpr
    uint64_extra(uint64_t const& a_v,
                 uint64_t const& a_extra)
        : v(a_v)
        , extra(a_extra)
    {}

    uint64_t v;
    uint64_t extra;
};

struct uint128_extra
{
    uint128_extra() = default;
    constexpr uint128_extra(uint128 const& a_v, uint64_t const& a_extra)
        : v(a_v)
        , extra(a_extra)
    {}

    uint128 v;
    uint64_t extra;
};

#endif

}  // namespace internals
}  // namespace softfloat

#endif  /* SOFTFLOAT_PRIMITIVES_TYPES_HPP_ */
