
/** @file

This C header file is part of the SoftFloat IEEE Floating-Point Arithmetic
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

#ifndef PRIMITIVES_FUNCTIONS_HPP_
#define PRIMITIVES_FUNCTIONS_HPP_

#include "primitives/types.hpp"

namespace softfloat {
namespace internals {

/**
Shifts `a' right by the number of bits given in `dist', which must be in
the range 1 to 63.  If any nonzero bits are shifted off, they are "jammed"
into the least-significant bit of the shifted value by setting the least-
significant bit to 1.  This shifted-and-jammed value is returned.
*/
inline constexpr uint64_t
short_shift_right_jam_64(uint64_t a,
                         uint8_t dist)
{
    return
        (a >> (dist >= 63 ? 63 : dist)) |
        ((((a << (63 - (dist >= 63 ? 63 : dist))) >> 1) + (~UINT64_C(0) >> 1)) >> 63);
}

/**
Shifts `a' right by the number of bits given in `dist', which must not
be zero.  If any nonzero bits are shifted off, they are "jammed" into the
least-significant bit of the shifted value by setting the least-significant
bit to 1.  This shifted-and-jammed value is returned.
  The value of `dist' can be arbitrarily large.  In particular, if `dist' is
greater than 32, the result will be either 0 or 1, depending on whether `a'
is zero or nonzero.
*/
inline uint32_t
shift_right_jam_32(uint32_t a,
                   uint16_t dist)
{
    auto const dist1 = dist >= 31 ? 31 : dist;
    auto const mask = ~(~UINT32_C(0) << dist1);
    auto const jam_bits = (a & mask) + mask;
    return (a | jam_bits) >> dist1;
}

/**
Shifts `a' right by the number of bits given in `dist', which must not
be zero.  If any nonzero bits are shifted off, they are "jammed" into the
least-significant bit of the shifted value by setting the least-significant
bit to 1.  This shifted-and-jammed value is returned.
  The value of `dist' can be arbitrarily large.  In particular, if `dist' is
greater than 64, the result will be either 0 or 1, depending on whether `a'
is zero or nonzero.
*/
inline uint64_t
shift_right_jam_64(uint64_t a,
                   uint32_t dist)
{
    auto const dist1 = dist >= 63 ? 63 : dist;
    auto const mask = ~(~UINT64_C(0) << dist1);
    auto const jam_bits = (a & mask) + mask;
    return (a | jam_bits) >> dist1;
}

/**
A constant table that translates an 8-bit unsigned integer (the array index)
into the number of leading 0 bits before the most-significant 1 of that
integer.  For integer zero (index 0), the corresponding table element is 8.
*/
extern const uint_least8_t count_leading_zeros_8[256];

/**
@returns the number of leading 0 bits before the most-significant 1 bit of `a'.  If `a' is zero, sizeof(Ty) * CHAR_BIT is returned.
*/
inline uint8_t
count_leading_zeros(uint8_t a)
{
    return count_leading_zeros_8[a];
}

inline uint8_t
count_leading_zeros(uint16_t a)
{
    return 0 == (a >> CHAR_BIT) ?
           count_leading_zeros(uint8_t(a)) + CHAR_BIT :
           count_leading_zeros(uint8_t(a >> CHAR_BIT));
}

inline uint8_t
count_leading_zeros(uint32_t a)
{
    return 0 == (a >> sizeof(uint16_t) * CHAR_BIT) ?
           count_leading_zeros(uint16_t(a)) + sizeof(uint16_t) * CHAR_BIT :
           count_leading_zeros(uint16_t(a >> sizeof(uint16_t) * CHAR_BIT));
}

inline uint8_t
count_leading_zeros(uint64_t a)
{
    return 0 == (a >> sizeof(uint32_t) * CHAR_BIT) ?
           count_leading_zeros(uint32_t(a)) + sizeof(uint32_t) * CHAR_BIT :
           count_leading_zeros(uint32_t(a >> sizeof(uint32_t) * CHAR_BIT));
}

extern const uint16_t approx_recip_1k0s[16];
extern const uint16_t approx_recip_1k1s[16];

/**
Returns an approximation to the reciprocal of the number represented by `a',
where `a' is interpreted as an unsigned fixed-point number with one integer
bit and 31 fraction bits.  The `a' input must be "normalized", meaning that
its most-significant bit (bit 31) must be 1.  Thus, if A is the value of
the fixed-point interpretation of `a', then 1 <= A < 2.  The returned value
is interpreted as a pure unsigned fraction, having no integer bits and 32
fraction bits.  The approximation returned is never greater than the true
reciprocal 1/A, and it differs from the true reciprocal by at most 2.006 ulp
(units in the last place).
*/

#if (SOFTFLOAT_FAST_DIV64TO32)

inline uint32_t
approx_recip_32_1(uint32_t a)
{
    return static_cast<uint32_t>(INT64_MAX / static_cast<uint32_t>(a));
}

#else

inline uint32_t
approx_recip_32_1(uint32_t a)
{
    auto const index = a >> 27 & 0xF;
    uint16_t const eps = static_cast<uint16_t>(a >> 11);
    uint16_t const r0 =
        approx_recip_1k0s[index] -
        ((approx_recip_1k1s[index] * static_cast<uint32_t>(eps)) >> 20);
    uint32_t const sigma0 = ~static_cast<uint32_t>((r0 * static_cast<uint64_t>(a)) >> 7);
    uint32_t const r =
        (static_cast<uint32_t>(r0) << 16) +
        static_cast<uint32_t>(static_cast<uint64_t>(r0) * sigma0 >> 24);
    uint32_t const sqrSigma0 = (static_cast<uint64_t>(sigma0) * sigma0) >> 32;
    return r + ((static_cast<uint32_t>(r) * static_cast<uint64_t>(sqrSigma0)) >> 48);
}

#endif

extern const uint16_t approx_recip_sqrt_1k0s[16];
extern const uint16_t approx_recip_sqrt_1k1s[16];

/**
Returns an approximation to the reciprocal of the square root of the number
represented by `a', where `a' is interpreted as an unsigned fixed-point
number either with one integer bit and 31 fraction bits or with two integer
bits and 30 fraction bits.  The format of `a' is determined by `oddExpA',
which must be either 0 or 1.  If `oddExpA' is 1, `a' is interpreted as
having one integer bit, and if `oddExpA' is 0, `a' is interpreted as having
two integer bits.  The `a' input must be "normalized", meaning that its
most-significant bit (bit 31) must be 1.  Thus, if A is the value of the
fixed-point interpretation of `a', it follows that 1 <= A < 2 when `oddExpA'
is 1, and 2 <= A < 4 when `oddExpA' is 0.
  The returned value is interpreted as a pure unsigned fraction, having
no integer bits and 32 fraction bits.  The approximation returned is never
greater than the true reciprocal 1/sqrt(A), and it differs from the true
reciprocal by at most 2.06 ulp (units in the last place).  The approximation
returned is also always within the range 0.5 to 1; thus, the most-
significant bit of the result is always set.
*/
uint32_t
approx_recip_sqrt_32_1(uint32_t oddExpA,
                       uint32_t a);

namespace fast_int64 {
using namespace internals;

/**
Returns true if the 128-bit unsigned integer formed by concatenating `a64'
and `a0' is less than or equal to the 128-bit unsigned integer formed by
concatenating `b64' and `b0'.
*/
inline constexpr bool
le(uint64_t const& a64,
   uint64_t const& a0,
   uint64_t const& b64,
   uint64_t const& b0)
{
    return a64 < b64 || (a64 == b64 && a0 <= b0);
}

inline constexpr bool
le(uint128 const& a,
   uint128 const& b)
{
    return le(a.v64, a.v0, b.v64, b.v64);
}

/**
Returns true if the 128-bit unsigned integer formed by concatenating `a64'
and `a0' is less than the 128-bit unsigned integer formed by concatenating
`b64' and `b0'.
*/
inline constexpr bool
lt(uint64_t const& a64,
   uint64_t const& a0,
   uint64_t const& b64,
   uint64_t const& b0)
{
    return
        a64 < b64 ||
        (a64 == b64 && a0 < b0);
}

inline constexpr bool
lt(uint128 const& a,
   uint128 const& b)
{
    return lt(a.v64, a.v0, b.v64, b.v0);
}

inline constexpr bool
lt(extFloat80_t const& a,
   extFloat80_t const& b)
{
    return lt(a.signExp, a.signif, b.signExp, b.signif);
}

/**
Shifts the 128 bits formed by concatenating `a64' and `a0' left by the
number of bits given in `dist', which must be in the range 1 to 63.
*/
inline constexpr uint128
short_shift_left_128(uint128 const& a,
                     uint8_t const dist)
{
    return uint128{a.v64 << dist | a.v0 >> (63u & -static_cast<int8_t>(dist)), a.v0 << dist};
}

/**
Shifts the 128 bits formed by concatenating `a64' and `a0' right by the
number of bits given in `dist', which must be in the range 1 to 63.
*/
inline constexpr uint128
short_shift_right_128(uint128 const& a,
                      uint8_t const dist)
{
    return uint128{a.v64 >> dist, a.v64 << (63 & -dist) | a.v0 >> dist};
}

/**
This function is the same as `shift_right_jam_64Extra' (below),
except that `dist' must be in the range 1 to 63.
*/
inline uint64_extra
shortShiftRightJam64Extra(uint64_t const& a,
                          uint64_t const& extra,
                          uint8_t const dist)
{
    return uint64_extra{a >> dist, a << (-static_cast<int8_t>(dist) & 63) | !!(0 != extra)};
}

/**
Shifts the 128 bits formed by concatenating `a64' and `a0' right by the
number of bits given in `dist', which must be in the range 1 to 63.  If any
nonzero bits are shifted off, they are "jammed" into the least-significant
bit of the shifted value by setting the least-significant bit to 1.  This
shifted-and-jammed value is returned.
*/
inline uint128
short_shift_right_jam_128(uint128 const& a,
                          uint8_t const dist)
{
    auto const negDist = 63 & -static_cast<int8_t>(dist);
    return uint128{
        a.v64 >> dist,
              a.v64 << negDist | a.v0 >> dist | !!(0 != (a.v0 << negDist))
    };
}

/**
This function is the same as `shift_right_jam_128Extra' (below),
except that `dist' must be in the range 1 to 63.
*/
inline uint128_extra
short_shift_right_jam_128Extra(uint128 const& a,
                               uint64_t const& extra,
                               uint8_t const dist)
{
    auto const uNegDist = 63 & -static_cast<int8_t>(dist);
    uint128_extra z;
    z.v.v64 = a.v64 >> dist;
    z.v.v0 = a.v64 << uNegDist | a.v0 >> dist;
    z.extra = a.v0 << uNegDist | (0 != extra);
    return z;
}

/**
Shifts the 128 bits formed by concatenating `a' and `extra' right by 64
_plus_ the number of bits given in `dist', which must not be zero.  This
shifted value is at most 64 nonzero bits and is returned in the `v' field
of the `struct uint64_extra' result.  The 64-bit `extra' field of the result
contains a value formed as follows from the bits that were shifted off:  The
_last_ bit shifted off is the most-significant bit of the `extra' field, and
the other 63 bits of the `extra' field are all zero if and only if _all_but_
_the_last_ bits shifted off were all zero.
(This function makes more sense if `a' and `extra' are considered to form
an unsigned fixed-point number with binary point between `a' and `extra'.
This fixed-point value is shifted right by the number of bits given in
`dist', and the integer part of this shifted value is returned in the `v'
field of the result.  The fractional part of the shifted value is modified
as described above and returned in the `extra' field of the result.)
*/
inline uint64_extra
shift_right_jam_64Extra(uint64_t a,
                        uint64_t extra,
                        uint32_t dist)
{
    uint64_extra z;

    if (dist < 64) {
        z.v = a >> dist;
        z.extra = a << (63 & -static_cast<int32_t>(dist));
    } else {
        z.v = 0;
        z.extra = dist == 64 ? a : !!(0 != a);
    }

    z.extra |= !!(extra != 0);
    return z;
}

/**
Shifts the 128 bits formed by concatenating `a64' and `a0' right by the
number of bits given in `dist', which must not be zero.  If any nonzero bits
are shifted off, they are "jammed" into the least-significant bit of the
shifted value by setting the least-significant bit to 1.  This shifted-and-
jammed value is returned.
The value of `dist' can be arbitrarily large.  In particular, if `dist' is
greater than 128, the result will be either 0 or 1, depending on whether the
original 128 bits are all zeros.
*/
uint128
shift_right_jam_128(uint64_t a64,
                    uint64_t a0,
                    uint32_t dist);

inline uint128
shift_right_jam_128(uint128 const& a,
                    uint32_t const dist)
{
    return shift_right_jam_128(a.v64, a.v0, dist);
}

/**
Shifts the 192 bits formed by concatenating `a64', `a0', and `extra' right
by 64 _plus_ the number of bits given in `dist', which must not be zero.
This shifted value is at most 128 nonzero bits and is returned in the `v'
field of the `struct uint128_extra' result.  The 64-bit `extra' field of the
result contains a value formed as follows from the bits that were shifted
off:  The _last_ bit shifted off is the most-significant bit of the `extra'
field, and the other 63 bits of the `extra' field are all zero if and only
if _all_but_the_last_ bits shifted off were all zero.
(This function makes more sense if `a64', `a0', and `extra' are considered
to form an unsigned fixed-point number with binary point between `a0' and
`extra'.  This fixed-point value is shifted right by the number of bits
given in `dist', and the integer part of this shifted value is returned
in the `v' field of the result.  The fractional part of the shifted value
is modified as described above and returned in the `extra' field of the
result.)
*/
uint128_extra
shift_right_jam_128Extra(uint64_t a64,
                         uint64_t a0,
                         uint64_t extra,
                         uint32_t dist);

inline uint128_extra
shift_right_jam_128Extra(uint128 const& a,
                         uint64_t const& extra,
                         uint32_t const& dist)
{
    return shift_right_jam_128Extra(a.v64, a.v0, extra, dist);
}

inline uint128_extra
shift_right_jam_128Extra(uint128_extra const& a,
                         uint32_t const& dist)
{
    return shift_right_jam_128Extra(a.v, a.extra, dist);
}

/**
Returns the sum of the 128-bit integer formed by concatenating `a64' and
`a0' and the 128-bit integer formed by concatenating `b64' and `b0'.  The
addition is modulo 2^128, so any carry out is lost.
*/
inline constexpr uint128
add(uint128 const& a,
    uint128 const& b)
{
    return uint128{a.v64 + b.v64 + !!(a.v0 + b.v0 < a.v0), a.v0 + b.v0};
}

/**
Adds the two 256-bit integers pointed to by `aPtr' and `bPtr'.  The addition
is modulo 2^256, so any carry out is lost.  The sum is stored at the
location pointed to by `zPtr'.  Each of `aPtr', `bPtr', and `zPtr' points to
an array of four 64-bit elements that concatenate in the platform's normal
endian order to form a 256-bit integer.
*/
inline void
add_M_256(uint64_t const* const aPtr,
          uint64_t const* const bPtr,
          uint64_t* const zPtr)
{
    bool carry = false;

    for (auto index = index_word_lo(4);; index += wordIncr) {
        uint64_t const wordA = aPtr[index];
        uint64_t const wordZ = wordA + bPtr[index] + !!carry;
        zPtr[index] = wordZ;

        if (index == index_word_hi(4)) {
            break;
        }

        if (wordZ != wordA) {
            carry = wordZ < wordA;
        }
    }
}

/**
Returns the difference of the 128-bit integer formed by concatenating `a64'
and `a0' and the 128-bit integer formed by concatenating `b64' and `b0'.
The subtraction is modulo 2^128, so any borrow out (carry out) is lost.
*/
inline constexpr uint128
sub(uint128 const& a,
    uint128 const& b)
{
    return uint128{a.v64 - b.v64 - !!(a.v0 < b.v0), a.v0 - b.v0};
}

/**
Subtracts the 256-bit integer pointed to by `bPtr' from the 256-bit integer
pointed to by `aPtr'.  The addition is modulo 2^256, so any borrow out
(carry out) is lost.  The difference is stored at the location pointed to
by `zPtr'.  Each of `aPtr', `bPtr', and `zPtr' points to an array of four
64-bit elements that concatenate in the platform's normal endian order to
form a 256-bit integer.
*/
void
sub_M_256(uint64_t const* aPtr,
          uint64_t const* bPtr,
          uint64_t* zPtr);

/**
@return the 128-bit product of `a', `b', and 2^32.
*/
inline uint128
mul_64_by_shifted_32_to_128(uint64_t const& a,
                            uint32_t const& b)
{
    uint64_t const mid = static_cast<uint64_t>(static_cast<uint32_t>(a)) * b;
    return
    uint128{
        mid << 32,
            static_cast<uint64_t>(static_cast<uint32_t>(a >> 32))* b + (mid >> 32)
    };
}

/**
@returns the 128-bit product of `a' and `b'.
*/
uint128
mul_64_to_128(uint64_t a,
              uint64_t b);

/**
@returns the product of the 128-bit integer formed by concatenating `a64' and
`a0', multiplied by `b'.  The multiplication is modulo 2^128; any overflow
bits are discarded.
*/
inline uint128
mul_128_by_32(uint64_t const& a64,
              uint64_t const& a0,
              uint32_t const& b)
{
    uint64_t const mid = static_cast<uint64_t>(static_cast<uint32_t>(a0 >> 32)) * b;
    uint64_t const v0 = a0 * b;
    uint32_t const carry = static_cast<uint32_t>(static_cast<uint32_t>(v0 >> 32) - static_cast<uint32_t>(mid));
    return uint128{a64* b + static_cast<uint32_t>((mid + carry) >> 32), v0};
}

inline uint128
mul_128_by_32(uint128 const& a,
              uint32_t const& b)
{
    return mul_128_by_32(a.v64, a.v0, b);
}

/**
Multiplies the 128-bit unsigned integer formed by concatenating `a64' and
`a0' by the 128-bit unsigned integer formed by concatenating `b64' and
`b0'.  The 256-bit product is stored at the location pointed to by `zPtr'.
Argument `zPtr' points to an array of four 64-bit elements that concatenate
in the platform's normal endian order to form a 256-bit integer.
*/
void
mul_M_128_to_256(uint64_t const& a64,
                 uint64_t const& a0,
                 uint64_t const& b64,
                 uint64_t const& b0,
                 uint64_t* zPtr);

inline void
mul_M_128_to_256(uint128 const& a,
                 uint128 const& b,
                 uint64_t* const zPtr)
{
    mul_M_128_to_256(a.v64, a.v0, b.v64, b.v0, zPtr);
}

}  // namespace fast_int64

namespace slow_int64 {

using namespace internals;

/**
The following functions are needed only when `SOFTFLOAT_FAST_INT64' is not
defined.
*/

/**
Compares the two 96-bit unsigned integers pointed to by `aPtr' and `bPtr'.
Returns -1 if the first integer (A) is less than the second (B); returns 0
if the two integers are equal; and returns +1 if the first integer (A)
is greater than the second (B).  (The result is thus the signum of A - B.)
Each of `aPtr' and `bPtr' points to an array of three 32-bit elements that
concatenate in the platform's normal endian order to form a 96-bit integer.
*/
int
compare_M_96(uint32_t const* aPtr,
             uint32_t const* bPtr);

/**
Compares the two 128-bit unsigned integers pointed to by `aPtr' and `bPtr'.
Returns -1 if the first integer (A) is less than the second (B); returns 0
if the two integers are equal; and returns +1 if the first integer (A)
is greater than the second (B).  (The result is thus the signum of A - B.)
Each of `aPtr' and `bPtr' points to an array of four 32-bit elements that
concatenate in the platform's normal endian order to form a 128-bit integer.
*/
int
compare_M_128(uint32_t const* aPtr,
              uint32_t const* bPtr);

/**
Extends `a' to 96 bits and shifts the value left by the number of bits given
in `dist', which must be in the range 1 to 31.  The result is stored at the
location pointed to by `zPtr'.  Argument `zPtr' points to an array of three
32-bit elements that concatenate in the platform's normal endian order to
form a 96-bit integer.
*/
inline void
short_shift_left_M_64_to_96(uint64_t a,
                           uint8_t dist,
                           uint32_t* zPtr)
{
    zPtr[index_word(3, 0)] = static_cast<uint32_t>(a) << dist;
    a >>= 32 - dist;
    zPtr[index_word(3, 2)] = a >> 32;
    zPtr[index_word(3, 1)] = static_cast<uint32_t>(a);
}

/**
Shifts the N-bit unsigned integer pointed to by `aPtr' left by the number
of bits given in `dist', where N = `size_words' * 32.  The value of `dist'
must be in the range 1 to 31.  Any nonzero bits shifted off are lost.  The
shifted N-bit result is stored at the location pointed to by `zPtr'.  Each
of `aPtr' and `zPtr' points to a `size_words'-long array of 32-bit elements
that concatenate in the platform's normal endian order to form an N-bit
integer.
*/
void
short_shift_left_M(size_t size_words,
                   uint32_t const* aPtr,
                   uint8_t dist,
                   uint32_t* zPtr);

/**
This function or macro is the same as `short_shift_left_M' with
`size_words' = 3 (N = 96).
*/
inline void
short_shift_left_M_96(uint32_t const* aPtr,
                     uint8_t dist,
                     uint32_t* zPtr)
{
    static constexpr size_t num_bits = 96u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    short_shift_left_M(num_words, aPtr, dist, zPtr);
}

/**
This function or macro is the same as `short_shift_left_M' with
`size_words' = 4 (N = 128).
*/
inline void
short_shift_left_M_128(uint32_t const* aPtr,
                       uint8_t dist,
                       uint32_t* zPtr)
{
    static constexpr size_t num_bits = 128u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    short_shift_left_M(num_words, aPtr, dist, zPtr);
}

/**
This function or macro is the same as `short_shift_left_M' with
`size_words' = 5 (N = 160).
*/
inline void
short_shift_left_M_160(uint32_t const* aPtr,
                       uint8_t dist,
                       uint32_t* zPtr)
{
    static constexpr size_t num_bits = 160u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    short_shift_left_M(num_words, aPtr, dist, zPtr);
}

/**
Shifts the N-bit unsigned integer left by the number of bits given in `dist', where N = `size_words' * 32.
Any nonzero bits shifted off are lost.
Each of `aPtr' and `zPtr' points to a `size_words'-long array of 32-bit elements that
concatenate in the platform's normal endian order to form an N-bit integer.

@param[in] size_words size of array in 32-bit elements
@param[in] aPtr pointer to N-bit unsigned integer
@param[in] dist shift in bits.
    The value of `dist' can be arbitrarily large.
    In particular, if `dist' is greater than N, the stored result will be 0.
    @pre The value of `dist' must not be zero.
@param[out] zPtr pointer to store location for the shifted N-bit result.
*/
void
shift_left_M(size_t size_words,
             uint32_t const* aPtr,
             uint32_t dist,
             uint32_t* zPtr);

/**
This function or macro is the same as `shift_left_M' with `size_words' = 3 (N = 96).
*/
inline void
shift_left_M_96(uint32_t const* aPtr,
                uint32_t dist,
                uint32_t* zPtr)
{
    static constexpr size_t num_bits = 96u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    shift_left_M(num_words, aPtr, dist, zPtr);
}

/**
This function or macro is the same as `shift_left_M' with `size_words' = 4 (N = 128).
*/
inline void
shift_left_M_128(uint32_t const* aPtr,
                 uint32_t dist,
                 uint32_t* zPtr)
{
    static constexpr size_t num_bits = 128u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    shift_left_M(num_words, aPtr, dist, zPtr);
}

/**
This function or macro is the same as `shift_left_M' with `size_words' = 5 (N = 160).
*/
inline void
shift_left_M_160(uint32_t const* aPtr,
                 uint32_t dist,
                 uint32_t* zPtr)
{
    static constexpr size_t num_bits = 160u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    shift_left_M(num_words, aPtr, dist, zPtr);
}

/**
Shifts the N-bit unsigned integer pointed to by `aPtr' right by the number
of bits given in `dist', where N = `size_words' * 32.  The value of `dist'
must be in the range 1 to 31.  Any nonzero bits shifted off are lost.  The
shifted N-bit result is stored at the location pointed to by `zPtr'.  Each
of `aPtr' and `zPtr' points to a `size_words'-long array of 32-bit elements
that concatenate in the platform's normal endian order to form an N-bit
integer.
*/
void
short_shift_right_M(size_t size_words,
                    uint32_t const* aPtr,
                    uint8_t dist,
                    uint32_t* zPtr);

/**
This function or macro is the same as `short_shift_right_M' with
`size_words' = 4 (N = 128).
*/
inline void
short_shift_right_M_128(uint32_t const* aPtr,
                        uint8_t dist,
                        uint32_t* zPtr)
{
    static constexpr size_t num_bits = 128u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    short_shift_right_M(num_words, aPtr, dist, zPtr);
}

/**
This function or macro is the same as `short_shift_right_M' with
`size_words' = 5 (N = 160).
*/
inline void
short_shift_right_M_160(uint32_t const* aPtr,
                        uint8_t dist,
                        uint32_t* zPtr)
{
    static constexpr size_t num_bits = 160u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    short_shift_right_M(num_words, aPtr, dist, zPtr);
}

/**
Shifts the N-bit unsigned integer pointed to by `aPtr' right by the number
of bits given in `dist', where N = `size_words' * 32.  The value of `dist'
must be in the range 1 to 31.  If any nonzero bits are shifted off, they are
"jammed" into the least-significant bit of the shifted value by setting the
least-significant bit to 1.  This shifted-and-jammed N-bit result is stored
at the location pointed to by `zPtr'.  Each of `aPtr' and `zPtr' points
to a `size_words'-long array of 32-bit elements that concatenate in the
platform's normal endian order to form an N-bit integer.
*/
void
short_shift_right_jam_M(size_t size_words,
                        uint32_t const* aPtr,
                        uint8_t dist,
                        uint32_t* zPtr);

/**
This function or macro is the same as `short_shift_right_jam_M' with
`size_words' = 5 (N = 160).
*/
inline void
short_shift_right_jam_M_160(uint32_t const* aPtr,
                            uint8_t dist,
                            uint32_t* zPtr)
{
    static constexpr size_t num_bits = 160u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    short_shift_right_jam_M(num_words, aPtr, dist, zPtr);
}

/**
Shifts the N-bit unsigned integer pointed to by `aPtr' right by the number
of bits given in `dist', where N = `size_words' * 32.  The value of `dist'
must not be zero.  Any nonzero bits shifted off are lost.  The shifted
N-bit result is stored at the location pointed to by `zPtr'.  Each of `aPtr'
and `zPtr' points to a `size_words'-long array of 32-bit elements that
concatenate in the platform's normal endian order to form an N-bit integer.
The value of `dist' can be arbitrarily large.  In particular, if `dist' is
greater than N, the stored result will be 0.
*/
void
shift_right_M(size_t size_words,
              uint32_t const* aPtr,
              uint32_t dist,
              uint32_t* zPtr);

/**
This function or macro is the same as `shift_right_M' with
`size_words' = 3 (N = 96).
*/
inline void
shift_right_M_96(uint32_t const* aPtr,
                 uint8_t dist,
                 uint32_t* zPtr)
{
    static constexpr size_t num_bits = 96u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    shift_right_M(num_words, aPtr, dist, zPtr);
}

/**
Shifts the N-bit unsigned integer pointed to by `aPtr' right by the number
of bits given in `dist', where N = `size_words' * 32.  The value of `dist'
must not be zero.  If any nonzero bits are shifted off, they are "jammed"
into the least-significant bit of the shifted value by setting the least-
significant bit to 1.  This shifted-and-jammed N-bit result is stored
at the location pointed to by `zPtr'.  Each of `aPtr' and `zPtr' points
to a `size_words'-long array of 32-bit elements that concatenate in the
platform's normal endian order to form an N-bit integer.
The value of `dist' can be arbitrarily large.  In particular, if `dist'
is greater than N, the stored result will be either 0 or 1, depending on
whether the original N bits are all zeros.
*/
void
shift_right_jam_M(size_t size_words,
                  uint32_t const* aPtr,
                  uint32_t dist,
                  uint32_t* zPtr);

/**
This function or macro is the same as `shift_right_jam_M' with
`size_words' = 3 (N = 96).
*/
inline void
shift_right_jam_M_96(uint32_t const* aPtr,
                     uint8_t dist,
                     uint32_t* zPtr)
{
    static constexpr size_t num_bits = 96u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    shift_right_jam_M(num_words, aPtr, dist, zPtr);
}

/**
This function or macro is the same as `shift_right_jam_M' with
`size_words' = 4 (N = 128).
*/
inline void
shift_right_jam_M_128(uint32_t const* aPtr,
                      uint8_t dist,
                      uint32_t* zPtr)
{
    static constexpr size_t num_bits = 128u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    shift_right_jam_M(num_words, aPtr, dist, zPtr);
}

/**
This function or macro is the same as `shift_right_jam_M' with
`size_words' = 5 (N = 160).
*/
inline void
shift_right_jam_M_160(uint32_t const* aPtr,
                      uint8_t dist,
                      uint32_t* zPtr)
{
    static constexpr size_t num_bits = 160u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    shift_right_jam_M(num_words, aPtr, dist, zPtr);
}

/**
Adds the two N-bit integers pointed to by `aPtr' and `bPtr', where N =
`size_words' * 32.  The addition is modulo 2^N, so any carry out is lost.
The N-bit sum is stored at the location pointed to by `zPtr'.  Each of
`aPtr', `bPtr', and `zPtr' points to a `size_words'-long array of 32-bit
elements that concatenate in the platform's normal endian order to form an
N-bit integer.
*/
void
add_M(size_t size_words,
      uint32_t const* aPtr,
      uint32_t const* bPtr,
      uint32_t* zPtr);

/**
This function or macro is the same as `add_M' with `size_words'
= 3 (N = 96).
*/
inline void
add_M_96(uint32_t const* aPtr,
         uint32_t const* bPtr,
         uint32_t* zPtr)
{
    static constexpr size_t num_bits = 96u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    add_M(num_words, aPtr, bPtr, zPtr);
}

/**
This function or macro is the same as `add_M' with `size_words'
= 4 (N = 128).
*/
inline void
add_M_128(uint32_t const* aPtr,
          uint32_t const* bPtr,
          uint32_t* zPtr)
{
    static constexpr size_t num_bits = 128u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    add_M(num_words, aPtr, bPtr, zPtr);
}

/**
This function or macro is the same as `add_M' with `size_words'
= 5 (N = 160).
*/
inline void
add_M_160(uint32_t const* aPtr,
          uint32_t const* bPtr,
          uint32_t* zPtr)
{
    static constexpr size_t num_bits = 160u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    add_M(num_words, aPtr, bPtr, zPtr);
}

/**
Adds the two N-bit unsigned integers pointed to by `aPtr' and `bPtr', where
N = `size_words' * 32, plus `carry', which must be either 0 or 1.  The N-bit
sum (modulo 2^N) is stored at the location pointed to by `zPtr', and any
carry out is returned as the result.  Each of `aPtr', `bPtr', and `zPtr'
points to a `size_words'-long array of 32-bit elements that concatenate in
the platform's normal endian order to form an N-bit integer.
*/
bool
add_carry_M(size_t size_words,
            uint32_t const* aPtr,
            uint32_t const* bPtr,
            bool carry,
            uint32_t* zPtr);

/**
This function or macro is the same as `add_carry_M', except that
the value of the unsigned integer pointed to by `bPtr' is bit-wise completed
before the addition.
*/
bool
add_compl_carry_M(size_t size_words,
                  uint32_t const* aPtr,
                  uint32_t const* bPtr,
                  bool carry,
                  uint32_t* zPtr);

/**
This function or macro is the same as `add_compl_carry_M' with
`size_words' = 3 (N = 96).
*/
inline bool
add_compl_carry_M_96(uint32_t const* aPtr,
                     uint32_t const* bPtr,
                     bool carry,
                     uint32_t* zPtr)
{
    static constexpr size_t num_bits = 96u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    return add_compl_carry_M(num_words, aPtr, bPtr, carry, zPtr);
}

/**
Replaces the N-bit unsigned integer pointed to by `zPtr' by the
2s-complement of itself, where N = `size_words' * 32.  Argument `zPtr'
points to a `size_words'-long array of 32-bit elements that concatenate in
the platform's normal endian order to form an N-bit integer.
*/
void
neg_M_X(uint8_t size_words,
        uint32_t* zPtr);

/**
This function or macro is the same as `neg_M_X' with `size_words'
= 3 (N = 96).
*/
inline void
neg_M_X96(uint32_t* zPtr)
{
    static constexpr size_t num_bits = 96u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    neg_M_X(num_words, zPtr);
}

/**
This function or macro is the same as `neg_M_X' with `size_words'
= 4 (N = 128).
*/
inline void
neg_M_X128(uint32_t* zPtr)
{
    static constexpr size_t num_bits = 128u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    neg_M_X(num_words, zPtr);
}

/**
This function or macro is the same as `neg_M_X' with `size_words' = 5 (N = 160).
*/
inline void
neg_M_X160(uint32_t* zPtr)
{
    static constexpr size_t num_bits = 160u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    neg_M_X(num_words, zPtr);
}

/**
This function or macro is the same as `neg_M_X' with `size_words' = 8 (N = 256).
*/
inline void
neg_M_X256(uint32_t* zPtr)
{
    static constexpr size_t num_bits = 256u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    neg_M_X(num_words, zPtr);
}

/**
Subtracts 1 from the N-bit integer pointed to by `zPtr', where N = `size_words' * 32.
The subtraction is modulo 2^N, so any borrow out (carry
out) is lost.  Argument `zPtr' points to a `size_words'-long array of 32-bit
elements that concatenate in the platform's normal endian order to form an
N-bit integer.
*/
void
sub_1_M_X(size_t size_words,
          uint32_t* zPtr);

/**
This function or macro is the same as `sub_1_M_X' with `size_words' = 3 (N = 96).
*/
inline void
sub_1_M_X96(uint32_t* zPtr)
{
    static constexpr size_t num_bits = 96u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    sub_1_M_X(num_words, zPtr);
}

/**
This function or macro is the same as `sub_1_M_X' with `size_words' = 5 (N = 160).
*/
inline void
sub_1_M_X160(uint32_t* zPtr)
{
    static constexpr size_t num_bits = 160u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    sub_1_M_X(num_words, zPtr);
}

/**
Subtracts the two N-bit integers pointed to by `aPtr' and `bPtr', where N = `size_words' * 32.
The subtraction is modulo 2^N, so any borrow out (carry out) is lost.
The N-bit difference is stored at the location pointed to by `zPtr'.
Each of `aPtr', `bPtr', and `zPtr' points to a `size_words'-long
array of 32-bit elements that concatenate in the platform's normal endian
order to form an N-bit integer.
*/
void
sub_M(size_t size_words,
      uint32_t const* aPtr,
      uint32_t const* bPtr,
      uint32_t* zPtr);

/**
This function or macro is the same as `sub_M' with `size_words' = 3 (N = 96).
*/
inline void
sub_M_96(uint32_t const* aPtr,
         uint32_t const* bPtr,
         uint32_t* zPtr)
{
    static constexpr size_t num_bits = 96u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    sub_M(num_words, aPtr, bPtr, zPtr);
}

/**
This function or macro is the same as `sub_M' with `size_words' = 4 (N = 128).
*/
inline void
sub_M_128(uint32_t const* aPtr,
          uint32_t const* bPtr,
          uint32_t* zPtr)
{
    static constexpr size_t num_bits = 128u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    sub_M(num_words, aPtr, bPtr, zPtr);
}

/**
This function or macro is the same as `sub_M' with `size_words' = 5 (N = 160).
*/
inline void
sub_M_160(uint32_t const* aPtr,
          uint32_t const* bPtr,
          uint32_t* zPtr)
{
    static constexpr size_t num_bits = 160u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    sub_M(num_words, aPtr, bPtr, zPtr);
}

/**
Multiplies `a' and `b' and stores the 128-bit product at the location
pointed to by `zPtr'.  Argument `zPtr' points to an array of four 32-bit
elements that concatenate in the platform's normal endian order to form a
128-bit integer.
*/
void
mul_M_64_to_128(uint64_t a,
                uint64_t b,
                uint32_t* zPtr);

/**
Multiplies the two 128-bit unsigned integers pointed to by `aPtr' and
`bPtr', and stores the 256-bit product at the location pointed to by `zPtr'.
Each of `aPtr' and `bPtr' points to an array of four 32-bit elements that
concatenate in the platform's normal endian order to form a 128-bit integer.
Argument `zPtr' points to an array of eight 32-bit elements that concatenate
to form a 256-bit integer.
*/
void
mul_M_128_to_256(uint32_t const* aPtr,
                 uint32_t const* bPtr,
                 uint32_t* zPtr);

/**
Performs a "remainder reduction step" as follows:  Arguments `remPtr' and
`bPtr' both point to N-bit unsigned integers, where N = `size_words' * 32.
Defining R and B as the values of those integers, the expression (R<<`dist')
- B * q is computed modulo 2^N, and the N-bit result is stored at the
location pointed to by `zPtr'.  Each of `remPtr', `bPtr', and `zPtr' points
to a `size_words'-long array of 32-bit elements that concatenate in the
platform's normal endian order to form an N-bit integer.
*/
void
rem_step_by_32_M(uint8_t size_words,
                 uint32_t const* remPtr,
                 uint8_t dist,
                 uint32_t const* bPtr,
                 uint32_t q,
                 uint32_t* zPtr);

/**
This function or macro is the same as `rem_step_by_32_M' with
`size_words' = 3 (N = 96).
*/
inline void
rem_step_by_32_M_96(uint32_t const* remPtr,
                    uint8_t dist,
                    uint32_t const* bPtr,
                    uint32_t q,
                    uint32_t* zPtr)
{
    static constexpr size_t num_bits = 96u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    rem_step_by_32_M(num_words, remPtr, dist, bPtr, q, zPtr);
}

/**
This function or macro is the same as `rem_step_by_32_M' with
`size_words' = 4 (N = 128).
*/
inline void
rem_step_by_32_M_128(uint32_t const* remPtr,
                     uint8_t dist,
                     uint32_t const* bPtr,
                     uint32_t q,
                     uint32_t* zPtr)
{
    static constexpr size_t num_bits = 128u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    rem_step_by_32_M(num_words, remPtr, dist, bPtr, q, zPtr);
}

/**
This function or macro is the same as `rem_step_by_32_M' with
`size_words' = 5 (N = 160).
*/
inline void
rem_step_M_160_by_32(uint32_t const* remPtr,
                     uint8_t dist,
                     uint32_t const* bPtr,
                     uint32_t q,
                     uint32_t* zPtr)
{
    static constexpr size_t num_bits = 160u;
    static_assert(0 == num_bits % (CHAR_BIT * sizeof(uint32_t)), "bad number of bits");
    static constexpr size_t num_words = num_bits / (CHAR_BIT * sizeof(uint32_t));
    rem_step_by_32_M(num_words, remPtr, dist, bPtr, q, zPtr);
}

}  // namespace slow_int64

}  // namespace internals
}  // namespace softfloat

#endif  /* PRIMITIVES_FUNCTIONS_HPP_ */
