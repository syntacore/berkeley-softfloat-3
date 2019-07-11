
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

#ifndef RISCV_TARGET_HPP_
#define RISCV_TARGET_HPP_

#include "softfloat/functions.h"
#include "primitives/types.hpp"
#include "internals.hpp"

namespace softfloat {
namespace internals {
inline namespace riscv {

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
static int32_t const i32_fromPosOverflow = INT32_MAX;
static int32_t const i32_fromNegOverflow = INT32_MIN;
static int32_t const i32_fromNaN = INT32_MAX;

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
The bit pattern for a default generated 16-bit floating-point NaN.
*/
static uint16_t const defaultNaNF16UI = UINT16_C(0x7E00);

/**
The bit pattern for a default generated 32-bit floating-point NaN.
*/
static uint32_t const defaultNaNF32UI = UINT32_C(0x7FC00000);

/**
"Common NaN" structure, used to transfer NaN representations from one format
to another.
*/
struct commonNaN
{
    commonNaN() = default;

    constexpr
    commonNaN(bool a_sign,
              uint64_t const& a_v64,
              uint64_t const& a_v0)
        : sign(a_sign)
        , v64(a_v64)
        , v0(a_v0)
    {}

private:
    bool sign;
    uint64_t v64;
    uint64_t v0;
};

/**
Assuming `uiA' has the bit pattern of a 16-bit floating-point NaN, converts
this NaN to the common NaN form, and stores the resulting common NaN at the
location pointed to by `zPtr'.  If the NaN is a signaling NaN, the invalid
exception is raised.
*/
inline commonNaN
softfloat_f16UIToCommonNaN(uint16_t const& uiA)
{
    if (0 == (uiA & 0x0200)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return commonNaN{false, 0u, 0u};
}

/**
Converts the common NaN pointed to by `aPtr' into a 16-bit floating-point
NaN, and returns the bit pattern of this value as an unsigned integer.
*/
inline constexpr uint16_t
softfloat_commonNaNToF16UI(commonNaN const&)
{
    return defaultNaNF16UI;
}

/**
Interpreting `uiA' and `uiB' as the bit patterns of two 16-bit floating-
point values, at least one of which is a NaN, returns the bit pattern of
the combined NaN result.  If either `uiA' or `uiB' has the pattern of a
signaling NaN, the invalid exception is raised.
*/
inline uint16_t
softfloat_propagateNaNF16UI(uint16_t const& uiA,
                            uint16_t const& uiB)
{
    if (softfloat_isSigNaNF16UI(uiA) || softfloat_isSigNaNF16UI(uiB)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return defaultNaNF16UI;
}

/**
Assuming `uiA' has the bit pattern of a 32-bit floating-point NaN, converts
this NaN to the common NaN form, and stores the resulting common NaN at the
location pointed to by `zPtr'.  If the NaN is a signaling NaN, the invalid
exception is raised.
*/
inline commonNaN
softfloat_f32UIToCommonNaN(uint32_t const& uiA)
{
    if (0 == (uiA & 0x00400000)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return commonNaN{false, 0u, 0u};
}

/**
Converts the common NaN pointed to by `aPtr' into a 32-bit floating-point
NaN, and returns the bit pattern of this value as an unsigned integer.
*/
inline constexpr uint32_t
softfloat_commonNaNToF32UI(commonNaN const&)
{
    return defaultNaNF32UI;
}

/**
The bit pattern for a default generated 64-bit floating-point NaN.
*/
static uint64_t const defaultNaNF64UI = UINT64_C(0x7FF8000000000000);

/**
Returns true when 64-bit unsigned integer `uiA' has the bit pattern of a
64-bit floating-point signaling NaN.
Note:  This macro evaluates its argument more than once.
*/
inline constexpr bool
softfloat_isSigNaNF64UI(uint64_t const& uiA)
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
inline commonNaN
softfloat_f64UIToCommonNaN(uint64_t const& uiA)
{
    if (0 == (uiA & UINT64_C(0x0008000000000000))) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return commonNaN{false, 0, 0};
}

/**
Converts the common NaN pointed to by `aPtr' into a 64-bit floating-point
NaN, and returns the bit pattern of this value as an unsigned integer.
*/
inline constexpr uint64_t
softfloat_commonNaNToF64UI(commonNaN const&)
{
    return defaultNaNF64UI;
}

/**
The bit pattern for a default generated 80-bit extended floating-point NaN.
*/
static uint16_t const defaultNaNExtF80UI64 = 0x7FFFu;
static uint64_t const defaultNaNExtF80UI0 = UINT64_C(0xC000000000000000);

/**
Returns true when the 80-bit unsigned integer formed from concatenating
16-bit `uiA64' and 64-bit `uiA0' has the bit pattern of an 80-bit extended
floating-point signaling NaN.
Note:  This macro evaluates its arguments more than once.
*/
inline constexpr bool
softfloat_isSigNaNExtF80UI(uint16_t const& uiA64,
                           uint64_t const& uiA0)
{
    return
        UINT16_C(0x7FFF) == (UINT16_C(0x7FFF) & uiA64) &&
        0 == (UINT64_C(0x4000000000000000) & uiA0) &&
        0 != (UINT64_C(0x3FFFFFFFFFFFFFFF) & uiA0);
}

template<typename Ty>
inline Ty
propagate_NaN(Ty const& uiA,
              Ty const& uiB);


/**
Interpreting `uiA' and `uiB' as the bit patterns of two 32-bit floating-
point values, at least one of which is a NaN, returns the bit pattern of
the combined NaN result.  If either `uiA' or `uiB' has the pattern of a
signaling NaN, the invalid exception is raised.
*/
template<>
inline uint32_t
propagate_NaN<uint32_t>(uint32_t const& uiA,
                        uint32_t const& uiB)
{
    if (softfloat_isSigNaNF32UI(uiA) || softfloat_isSigNaNF32UI(uiB)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return defaultNaNF32UI;
}

/**
Interpreting `uiA' and `uiB' as the bit patterns of two 64-bit floating-
point values, at least one of which is a NaN, returns the bit pattern of
the combined NaN result.  If either `uiA' or `uiB' has the pattern of a
signaling NaN, the invalid exception is raised.
*/
static inline uint64_t
softfloat_propagateNaNF64UI(uint64_t const& uiA,
                            uint64_t const& uiB)
{
    if (softfloat_isSigNaNF64UI(uiA) || softfloat_isSigNaNF64UI(uiB)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return defaultNaNF64UI;
}

#ifdef SOFTFLOAT_FAST_INT64

/*
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
inline commonNaN
softfloat_extF80UIToCommonNaN(uint64_t const&,
                              uint64_t const& uiA0)
{
    if (0 == (uiA0 & UINT64_C(0x4000000000000000))) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return commonNaN{false, 0, 0};
}

/** The bit pattern for a default generated 128-bit floating-point NaN. */
/**@{*/
static uint64_t const defaultNaNF128UI64 = UINT64_C(0x7FFF800000000000);
static uint64_t const defaultNaNF128UI0 = UINT64_C(0);
/**@}*/

/**
Returns true when the 128-bit unsigned integer formed from concatenating
64-bit `uiA64' and 64-bit `uiA0' has the bit pattern of a 128-bit floating-
point signaling NaN.
Note:  This macro evaluates its arguments more than once.
*/
inline constexpr bool
softfloat_isSigNaNF128UI(uint64_t const& uiA64,
                         uint64_t const& uiA0)
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
inline commonNaN
softfloat_f128UIToCommonNaN(uint64_t const& uiA64,
                            uint64_t const&)
{
    if (0 == (uiA64 & UINT64_C(0x0000800000000000))) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return commonNaN{false, 0, 0};
}

/**
Converts the common NaN pointed to by `aPtr' into an 80-bit extended
floating-point NaN, and returns the bit pattern of this value as an unsigned
integer.
*/
inline constexpr uint128
softfloat_commonNaNToExtF80UI(commonNaN const&)
{
    return uint128{defaultNaNExtF80UI64, defaultNaNExtF80UI0};
}

/**
Interpreting the unsigned integer formed from concatenating `uiA64' and
`uiA0' as an 80-bit extended floating-point value, and likewise interpreting
the unsigned integer formed from concatenating `uiB64' and `uiB0' as another
80-bit extended floating-point value, and assuming at least on of these
floating-point values is a NaN, returns the bit pattern of the combined NaN
result.  If either original floating-point value is a signaling NaN, the
invalid exception is raised.
*/
inline uint128
softfloat_propagateNaNExtF80UI(uint16_t const& uiA64,
                               uint64_t const& uiA0,
                               uint16_t const& uiB64,
                               uint64_t const& uiB0)
{

    if (softfloat_isSigNaNExtF80UI(uiA64, uiA0) || softfloat_isSigNaNExtF80UI(uiB64, uiB0)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return uint128{defaultNaNExtF80UI64, defaultNaNExtF80UI0};
}

/**
Converts the common NaN pointed to by `aPtr' into a 128-bit floating-point
NaN, and returns the bit pattern of this value as an unsigned integer.
*/
inline constexpr uint128
softfloat_commonNaNToF128UI(commonNaN const&)
{
    return uint128{defaultNaNF128UI64, defaultNaNF128UI0};
}

/**
Interpreting the unsigned integer formed from concatenating `uiA64' and
`uiA0' as a 128-bit floating-point value, and likewise interpreting the
unsigned integer formed from concatenating `uiB64' and `uiB0' as another
128-bit floating-point value, and assuming at least on of these floating-
point values is a NaN, returns the bit pattern of the combined NaN result.
If either original floating-point value is a signaling NaN, the invalid
exception is raised.
*/
inline uint128
softfloat_propagateNaNF128UI(uint64_t const& uiA64,
                             uint64_t const& uiA0,
                             uint64_t const& uiB64,
                             uint64_t const& uiB0)
{
    if (softfloat_isSigNaNF128UI(uiA64, uiA0) || softfloat_isSigNaNF128UI(uiB64, uiB0)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return uint128{defaultNaNF128UI64, defaultNaNF128UI0};
}

#else

/**
The bit pattern for a default generated 128-bit floating-point NaN.
*/
static uint32_t const defaultNaNF128UI96 = UINT32_C(0x7FFF8000);
static uint32_t const defaultNaNF128UI64 = UINT32_C(0);
static uint32_t const defaultNaNF128UI32 = UINT32_C(0);
static uint32_t const defaultNaNF128UI0 = UINT32_C(0);

/**
Assuming the 80-bit extended floating-point value pointed to by `aSPtr' is
a NaN, converts this NaN to the common NaN form, and stores the resulting
common NaN at the location pointed to by `zPtr'.  If the NaN is a signaling
NaN, the invalid exception is raised.
*/
inline commonNaN
softfloat_extF80MToCommonNaN(extFloat80_t const& a)
{
    if (0 == (a.signif & UINT64_C(0x4000000000000000))) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return commonNaN{false, 0, 0};
}

/**
Assuming the 128-bit floating-point value pointed to by `aWPtr' is a NaN,
converts this NaN to the common NaN form, and stores the resulting common
NaN at the location pointed to by `zPtr'.  If the NaN is a signaling NaN,
the invalid exception is raised.  Argument `aWPtr' points to an array of
four 32-bit elements that concatenate in the platform's normal endian order
to form a 128-bit floating-point value.
*/
inline commonNaN
softfloat_f128MToCommonNaN(uint32_t const* aWPtr)
{
    if (0 == (aWPtr[indexWordHi(4)] & UINT64_C(0x0000800000000000))) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    return commonNaN{false, 0, 0};
}

/**
Converts the common NaN pointed to by `aPtr' into an 80-bit extended
floating-point NaN, and stores this NaN at the location pointed to by
`zSPtr'.
*/
/** @bug use extFloat80_t */
inline extFloat80M
softfloat_commonNaNToExtF80M(commonNaN const&)
{
    extFloat80M z;
    z.signExp = defaultNaNExtF80UI64;
    z.signif = defaultNaNExtF80UI0;
    return z;
}

/**
Assuming at least one of the two 80-bit extended floating-point values
pointed to by `aSPtr' and `bSPtr' is a NaN, stores the combined NaN result
at the location pointed to by `zSPtr'.  If either original floating-point
value is a signaling NaN, the invalid exception is raised.
*/
/** @bug use extFloat80_t */
inline void
softfloat_propagateNaNExtF80M(extFloat80M const* const aSPtr,
                              extFloat80M const* const bSPtr,
                              extFloat80M* const zSPtr)
{
    if (softfloat_isSigNaNExtF80UI(aSPtr->signExp, aSPtr->signif) || (bSPtr && softfloat_isSigNaNExtF80UI(bSPtr->signExp, bSPtr->signif))) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    zSPtr->signExp = defaultNaNExtF80UI64;
    zSPtr->signif = defaultNaNExtF80UI0;
}

/**
Converts the common NaN pointed to by `aPtr' into a 128-bit floating-point
NaN, and stores this NaN at the location pointed to by `zWPtr'.  Argument
`zWPtr' points to an array of four 32-bit elements that concatenate in the
platform's normal endian order to form a 128-bit floating-point value.
*/
inline void
softfloat_commonNaNToF128M(commonNaN const&,
                           uint32_t* const zWPtr)
{
    zWPtr[indexWord(4, 3)] = defaultNaNF128UI96;
    zWPtr[indexWord(4, 2)] = defaultNaNF128UI64;
    zWPtr[indexWord(4, 1)] = defaultNaNF128UI32;
    zWPtr[indexWord(4, 0)] = defaultNaNF128UI0;
}

/**
Assuming at least one of the two 128-bit floating-point values pointed to by
`aWPtr' and `bWPtr' is a NaN, stores the combined NaN result at the location
pointed to by `zWPtr'.  If either original floating-point value is a
signaling NaN, the invalid exception is raised.  Each of `aWPtr', `bWPtr',
and `zWPtr' points to an array of four 32-bit elements that concatenate in
the platform's normal endian order to form a 128-bit floating-point value.
*/
inline void
softfloat_propagateNaNF128M(uint32_t const* const aWPtr,
                            uint32_t const* const bWPtr,
                            uint32_t* const zWPtr)
{
    if (f128M_isSignalingNaN((const float128_t*)aWPtr) || (bWPtr && f128M_isSignalingNaN((const float128_t*)bWPtr))) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    zWPtr[indexWord(4, 3)] = defaultNaNF128UI96;
    zWPtr[indexWord(4, 2)] = defaultNaNF128UI64;
    zWPtr[indexWord(4, 1)] = defaultNaNF128UI32;
    zWPtr[indexWord(4, 0)] = defaultNaNF128UI0;
}

#endif

}  // namespace riscv
}  // namespace internals
}  // namespace softfloat

#endif  /* RISCV_TARGET_HPP_ */
