/** @file

This C header file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

@copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
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

#ifndef SPECIALIZE_H_
#define SPECIALIZE_H_

#include "softfloat/functions.h"
#include "primitives/types.hpp"
#include "internals.hpp"

namespace softfloat {
namespace internals {
inline namespace Intel_8086 {

/**
Default value for `softfloat_detectTininess'.
*/
static softfloat_tininess const init_detectTininess = softfloat_tininess_afterRounding;

/**
The values to return on conversions to 32-bit integer formats that raise an
invalid exception.
*/
static uint32_t const ui32_fromPosOverflow = UINT32_MAX;
static uint32_t const ui32_fromNegOverflow = 0;
static uint32_t const ui32_fromNaN = UINT32_MAX;
static int32_t const i32_fromPosOverflow  = INT32_MAX;
static int32_t const i32_fromNegOverflow  = INT32_MIN;
static int32_t const i32_fromNaN          = INT32_MAX;

/**
The values to return on conversions to 64-bit integer formats that raise an
invalid exception.
*/
static uint64_t const ui64_fromPosOverflow = UINT64_MAX;
static uint64_t const ui64_fromNegOverflow = 0;
static uint64_t const ui64_fromNaN = UINT64_MAX;
static int64_t const i64_fromPosOverflow = INT64_MAX;
static int64_t const i64_fromNegOverflow = INT64_MIN;
static int64_t const i64_fromNaN = INT64_MAX;

/**
"Common NaN" structure, used to transfer NaN representations from one format
to another.
*/
struct commonNaN {
    bool sign;
#ifdef BIG_ENDIAN
    uint64_t v64, v0;
#else
    uint64_t v0, v64;
#endif
};

/**
The bit pattern for a default generated 16-bit floating-point NaN.
*/
static uint16_t const defaultNaNF16UI = UINT16_C(0xFE00);

/**
Returns true when 16-bit unsigned integer `uiA' has the bit pattern of a
16-bit floating-point signaling NaN.
Note:  This macro evaluates its argument more than once.
*/
inline bool
softfloat_isSigNaNF16UI(uint16_t uiA)
{
    return UINT16_C(0x7C00) == (uiA & UINT16_C(0x7E00)) && 0 != (uiA & 0x01FF);
}

/**
Assuming `uiA' has the bit pattern of a 16-bit floating-point NaN, converts
this NaN to the common NaN form, and stores the resulting common NaN at the
location pointed to by `zPtr'.  If the NaN is a signaling NaN, the invalid
exception is raised.
*/
commonNaN
softfloat_f16UIToCommonNaN(uint16_t uiA);

/**
Converts the common NaN pointed to by `aPtr' into a 16-bit floating-point
NaN, and returns the bit pattern of this value as an unsigned integer.
*/
inline constexpr uint16_t
softfloat_commonNaNToF16UI(commonNaN a)
{
    return static_cast<uint16_t>(a.sign) << 15 | 0x7E00 | a.v64 >> 54;
}

/**
Interpreting `uiA' and `uiB' as the bit patterns of two 16-bit floating-
point values, at least one of which is a NaN, returns the bit pattern of
the combined NaN result.  If either `uiA' or `uiB' has the pattern of a
signaling NaN, the invalid exception is raised.
*/
uint16_t
softfloat_propagateNaNF16UI(uint16_t uiA, uint16_t uiB);

/**
The bit pattern for a default generated 32-bit floating-point NaN.
*/
uint32_t const defaultNaNF32UI = 0xFFC00000u;

/**
Assuming `uiA' has the bit pattern of a 32-bit floating-point NaN, converts
this NaN to the common NaN form, and stores the resulting common NaN at the
location pointed to by `zPtr'.  If the NaN is a signaling NaN, the invalid
exception is raised.
*/
commonNaN
softfloat_f32UIToCommonNaN(uint32_t uiA);

/**
Converts the common NaN pointed to by `aPtr' into a 32-bit floating-point
NaN, and returns the bit pattern of this value as an unsigned integer.
*/
inline constexpr uint32_t
softfloat_commonNaNToF32UI(commonNaN a)
{
    return static_cast<uint32_t>(a.sign) << 31 | 0x7FC00000 | a.v64 >> 41;
}

/**
Interpreting `uiA' and `uiB' as the bit patterns of two 32-bit floating-
point values, at least one of which is a NaN, returns the bit pattern of
the combined NaN result.  If either `uiA' or `uiB' has the pattern of a
signaling NaN, the invalid exception is raised.
*/
uint32_t
softfloat_propagateNaNF32UI(uint32_t uiA,
                            uint32_t uiB);

/**
The bit pattern for a default generated 64-bit floating-point NaN.
*/
uint64_t const defaultNaNF64UI = UINT64_C(0xFFF8000000000000);

/**
Returns true when 64-bit unsigned integer `uiA' has the bit pattern of a
64-bit floating-point signaling NaN.
Note:  This macro evaluates its argument more than once.
*/
inline constexpr bool
softfloat_isSigNaNF64UI(uint64_t uiA)
{
    return
        UINT64_C(0x7FF0000000000000) == (uiA & UINT64_C(0x7FF8000000000000)) &&
        0 != (uiA & UINT64_C(0x0007FFFFFFFFFFFF));
}

/**
Assuming `uiA' has the bit pattern of a 64-bit floating-point NaN, converts
this NaN to the common NaN form, and stores the resulting common NaN at the
location pointed to by `zPtr'.  If the NaN is a signaling NaN, the invalid
exception is raised.
*/
commonNaN
softfloat_f64UIToCommonNaN(uint64_t uiA);

/**
Converts the common NaN pointed to by `aPtr' into a 64-bit floating-point
NaN, and returns the bit pattern of this value as an unsigned integer.
*/
inline constexpr uint64_t
softfloat_commonNaNToF64UI(commonNaN a)
{
    return
        static_cast<uint64_t>(a.sign) << 63 | UINT64_C(0x7FF8000000000000) |
        a.v64 >> 12;
}

/**
Interpreting `uiA' and `uiB' as the bit patterns of two 64-bit floating-
point values, at least one of which is a NaN, returns the bit pattern of
the combined NaN result.  If either `uiA' or `uiB' has the pattern of a
signaling NaN, the invalid exception is raised.
*/
uint64_t
softfloat_propagateNaNF64UI(uint64_t uiA,
                            uint64_t uiB);

/**
The bit pattern for a default generated 80-bit extended floating-point NaN.
*/
uint16_t const defaultNaNExtF80UI64 = 0xFFFFu;
uint64_t const defaultNaNExtF80UI0 = UINT64_C(0xC000000000000000);

/**
Returns true when the 80-bit unsigned integer formed from concatenating
16-bit `uiA64' and 64-bit `uiA0' has the bit pattern of an 80-bit extended
floating-point signaling NaN.
Note:  This macro evaluates its arguments more than once.
*/
inline constexpr bool
softfloat_isSigNaNExtF80UI(uint16_t uiA64, uint64_t uiA0)
{
    return
        0x7FFF == (uiA64 & 0x7FFF) &&
        0 == (uiA0 & UINT64_C(0x4000000000000000)) &&
        0 != (uiA0 & UINT64_C(0x3FFFFFFFFFFFFFFF));
}

#ifdef SOFTFLOAT_FAST_INT64

/**
The following functions are needed only when `SOFTFLOAT_FAST_INT64' is
defined.
*/

/**
Assuming the unsigned integer formed from concatenating `uiA64' and `uiA0'
has the bit pattern of an 80-bit extended floating-point NaN, converts
this NaN to the common NaN form, and stores the resulting common NaN at the
location pointed to by `zPtr'.  If the NaN is a signaling NaN, the invalid
exception is raised.
*/
commonNaN
softfloat_extF80UIToCommonNaN(uint16_t uiA64,
                              uint64_t uiA0);

/**
Converts the common NaN pointed to by `aPtr' into an 80-bit extended
floating-point NaN, and returns the bit pattern of this value as an unsigned
integer.
*/
uint128
softfloat_commonNaNToExtF80UI(commonNaN a);

/**
Interpreting the unsigned integer formed from concatenating `uiA64' and
`uiA0' as an 80-bit extended floating-point value, and likewise interpreting
the unsigned integer formed from concatenating `uiB64' and `uiB0' as another
80-bit extended floating-point value, and assuming at least on of these
floating-point values is a NaN, returns the bit pattern of the combined NaN
result.  If either original floating-point value is a signaling NaN, the
invalid exception is raised.
*/
uint128
softfloat_propagateNaNExtF80UI(uint16_t uiA64,
                               uint64_t uiA0,
                               uint16_t uiB64,
                               uint64_t uiB0);

/**
The bit pattern for a default generated 128-bit floating-point NaN.
*/
static uint64_t const defaultNaNF128UI64 = UINT64_C(0xFFFF800000000000);
static uint64_t const defaultNaNF128UI0 = UINT64_C(0);

/**
Returns true when the 128-bit unsigned integer formed from concatenating
64-bit `uiA64' and 64-bit `uiA0' has the bit pattern of a 128-bit floating-
point signaling NaN.
Note:  This macro evaluates its arguments more than once.
*/
inline constexpr bool
softfloat_isSigNaNF128UI(uint64_t uiA64,
                         uint64_t uiA0)
{
    return
        UINT64_C(0x7FFF000000000000) == (uiA64 & UINT64_C(0x7FFF800000000000)) &&
        (0 != uiA0 || 0 != (uiA64 & UINT64_C(0x00007FFFFFFFFFFF)));
}

/**
Assuming the unsigned integer formed from concatenating `uiA64' and `uiA0'
has the bit pattern of a 128-bit floating-point NaN, converts this NaN to
the common NaN form, and stores the resulting common NaN at the location
pointed to by `zPtr'.  If the NaN is a signaling NaN, the invalid exception
is raised.
*/
commonNaN
softfloat_f128UIToCommonNaN(uint64_t uiA64,
                            uint64_t uiA0);

/**
Converts the common NaN pointed to by `aPtr' into a 128-bit floating-point
NaN, and returns the bit pattern of this value as an unsigned integer.
*/
uint128
softfloat_commonNaNToF128UI(commonNaN);

/**
Interpreting the unsigned integer formed from concatenating `uiA64' and
`uiA0' as a 128-bit floating-point value, and likewise interpreting the
unsigned integer formed from concatenating `uiB64' and `uiB0' as another
128-bit floating-point value, and assuming at least on of these floating-
point values is a NaN, returns the bit pattern of the combined NaN result.
If either original floating-point value is a signaling NaN, the invalid
exception is raised.
*/
uint128
softfloat_propagateNaNF128UI(uint64_t uiA64,
                             uint64_t uiA0,
                             uint64_t uiB64,
                             uint64_t uiB0);

#else

/**
The following functions are needed only when `SOFTFLOAT_FAST_INT64' is not
defined.
*/

/**
Assuming the 80-bit extended floating-point value pointed to by `aSPtr' is
a NaN, converts this NaN to the common NaN form, and stores the resulting
common NaN at the location pointed to by `zPtr'.  If the NaN is a signaling
NaN, the invalid exception is raised.
*/
/** @bug use extFloat80_t */
commonNaN
softfloat_extF80MToCommonNaN(extFloat80M aSPtr);

/**
Converts the common NaN pointed to by `aPtr' into an 80-bit extended
floating-point NaN, and stores this NaN at the location pointed to by
`zSPtr'.
*/
/** @bug use extFloat80_t */
extFloat80M
softfloat_commonNaNToExtF80M(commonNaN a);

/**
Assuming at least one of the two 80-bit extended floating-point values
pointed to by `aSPtr' and `bSPtr' is a NaN, stores the combined NaN result
at the location pointed to by `zSPtr'.  If either original floating-point
value is a signaling NaN, the invalid exception is raised.
*/
/** @bug use extFloat80_t */
void
softfloat_propagateNaNExtF80M(extFloat80M const* aSPtr,
                              extFloat80M const* bSPtr,
                              extFloat80M* zSPtr);

/**
The bit pattern for a default generated 128-bit floating-point NaN.
*/
static uint32_t const defaultNaNF128UI96 = UINT32_C(0xFFFF8000);
static uint32_t const defaultNaNF128UI64 = UINT32_C(0);
static uint32_t const defaultNaNF128UI32 = UINT32_C(0);
static uint32_t const defaultNaNF128UI0  = UINT32_C(0);

/**
Assuming the 128-bit floating-point value pointed to by `aWPtr' is a NaN,
converts this NaN to the common NaN form, and stores the resulting common
NaN at the location pointed to by `zPtr'.  If the NaN is a signaling NaN,
the invalid exception is raised.  Argument `aWPtr' points to an array of
four 32-bit elements that concatenate in the platform's normal endian order
to form a 128-bit floating-point value.
*/
commonNaN
softfloat_f128MToCommonNaN(uint32_t const* aWPtr);

/**
Converts the common NaN pointed to by `aPtr' into a 128-bit floating-point
NaN, and stores this NaN at the location pointed to by `zWPtr'.  Argument
`zWPtr' points to an array of four 32-bit elements that concatenate in the
platform's normal endian order to form a 128-bit floating-point value.
*/
void
softfloat_commonNaNToF128M(commonNaN a,
                           uint32_t* zWPtr);

/**
Assuming at least one of the two 128-bit floating-point values pointed to by
`aWPtr' and `bWPtr' is a NaN, stores the combined NaN result at the location
pointed to by `zWPtr'.  If either original floating-point value is a
signaling NaN, the invalid exception is raised.  Each of `aWPtr', `bWPtr',
and `zWPtr' points to an array of four 32-bit elements that concatenate in
the platform's normal endian order to form a 128-bit floating-point value.
*/
void
softfloat_propagateNaNF128M(uint32_t const* aWPtr,
                            uint32_t const* bWPtr,
                            uint32_t* zWPtr);

#endif  /* SOFTFLOAT_FAST_INT64 */
}  // namespace Intel_8086
}  // namespace internals
}  // namespace softfloat

#endif  /* SPECIALIZE_H_ */
