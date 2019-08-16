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

#ifndef INTEL8086_SSE_TARGET_HPP_
#define INTEL8086_SSE_TARGET_HPP_

#include "Intel8086.hpp"

namespace softfloat {
namespace internals {
inline namespace Intel_8086 {

/**
Interpreting `uiA' and `uiB' as the bit patterns of two 16-bit floating-
point values, at least one of which is a NaN, returns the bit pattern of
the combined NaN result.  If either `uiA' or `uiB' has the pattern of a
signaling NaN, the invalid exception is raised.
*/
inline uint16_t
softfloat_propagateNaNF16UI(uint16_t const uiA,
                            uint16_t const uiB)
{
    static uint16_t const quite_nan_bit = UINT16_C(0x0200);
    bool const isSigNaNA = is_sNaN(uiA);

    if (isSigNaNA || is_sNaN(uiB)) {
        softfloat_raiseFlags(softfloat_flag_invalid);

        if (isSigNaNA) {
            return static_cast<uint16_t>(quite_nan_bit | uiA);
        }
    }

    return static_cast<uint16_t>(quite_nan_bit | (is_NaN(uiA) ? uiA : uiB));
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
    static uint32_t const quietNaN_bit = UINT32_C(0x00400000);

    bool const isSigNaNA = softfloat_isSigNaNF32UI(uiA);
    bool const isSigNaNB = softfloat_isSigNaNF32UI(uiB);
    /* Make NaNs non-signaling. */
    uint32_t const uiNonsigA = quietNaN_bit | uiA;
    uint32_t const uiNonsigB = quietNaN_bit | uiB;

    if (isSigNaNA || isSigNaNB) {
        softfloat_raiseFlags(softfloat_flag_invalid);

        if (!isSigNaNA) {
            return is_NaN(uiA) ? uiNonsigA : uiNonsigB;
        }

        if (!isSigNaNB) {
            return is_NaN(uiB) ? uiNonsigB : uiNonsigA;
        }
    }

    {
        uint32_t const uiMagA = uiNonsigA & 0x7FFFFFFF;
        uint32_t const uiMagB = uiNonsigB & 0x7FFFFFFF;
        return
            uiMagA < uiMagB ? uiNonsigB :
            uiMagB < uiMagA ? uiNonsigA :
            uiNonsigA < uiNonsigB ? uiNonsigA : uiNonsigB;
    }
}

/**
Interpreting `uiA' and `uiB' as the bit patterns of two 64-bit floating-
point values, at least one of which is a NaN, returns the bit pattern of
the combined NaN result.  If either `uiA' or `uiB' has the pattern of a
signaling NaN, the invalid exception is raised.
*/
inline uint64_t
propagate_NaN(uint64_t const uiA,
                            uint64_t const uiB)
{
    bool const isSigNaNA = is_sNaN(uiA);

    if (isSigNaNA || is_sNaN(uiB)) {
        softfloat_raiseFlags(softfloat_flag_invalid);

        if (isSigNaNA) {
            return uiA | UINT64_C(0x0008000000000000);
        }
    }

    return (is_NaN(uiA) ? uiA : uiB) | UINT64_C(0x0008000000000000);
}

}  // namespace Intel_8086
}  // namespace internals
}  // namespace softfloat

#endif  /* INTEL8086_SSE_TARGET_HPP_ */
