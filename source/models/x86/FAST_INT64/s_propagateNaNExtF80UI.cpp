
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

#ifndef SOFTFLOAT_FAST_INT64
#error Fast int64_t operations only
#endif

namespace softfloat {
namespace internals {
namespace Intel_8086 {

namespace {
static uint128
return_larger_magnitude(uint16_t const uiA64,
                        uint64_t const uiNonsigA0,
                        uint16_t const uiB64,
                        uint64_t const uiNonsigB0)
{
    uint16_t const uiMagA64 = uiA64 & 0x7FFFu;
    uint16_t const uiMagB64 = uiB64 & 0x7FFFu;

    if (uiMagA64 < uiMagB64) {
        return uint128{uiB64, uiNonsigB0};
    }

    if (uiMagB64 < uiMagA64) {
        return uint128{uiA64, uiNonsigA0};
    }

    assert(uiMagA64 == uiMagB64);

    if (uiNonsigA0 < uiNonsigB0) {
        return uint128{uiB64, uiNonsigB0};
    }

    if (uiNonsigB0 < uiNonsigA0) {
        return uint128{uiA64, uiNonsigA0};
    }

    assert(uiNonsigA0 == uiNonsigB0);

    return
        uiA64 < uiB64 ?
        uint128{uiA64, uiNonsigA0} :
        uint128{uiB64, uiNonsigB0};
}

}  // namespace

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
softfloat_propagateNaNExtF80UI(uint16_t const uiA64,
                               uint64_t const uiA0,
                               uint16_t const uiB64,
                               uint64_t const uiB0)
{
    bool const isSigNaNA = is_sNaN(uiA64, uiA0);
    bool const isSigNaNB = is_sNaN(uiB64, uiB0);
    /* Make NaNs non-signaling. */
    uint64_t const uiNonsigA0 = uiA0 | UINT64_C(0xC000000000000000);
    uint64_t const uiNonsigB0 = uiB0 | UINT64_C(0xC000000000000000);

    if (isSigNaNA || isSigNaNB) {
        softfloat_raiseFlags(softfloat_flag_invalid);

        if (isSigNaNA) {
            if (isSigNaNB) {
                return return_larger_magnitude(uiA64, uiNonsigA0, uiB64, uiNonsigB0);
            }

            if (is_NaN(uiB64, uiB0)) {
                return uint128{uiB64, uiNonsigB0};
            }

            return uint128{uiA64, uiNonsigA0};
        }

        if (is_NaN(uiA64, uiA0)) {
            return uint128{uiA64, uiNonsigA0};
        }

        return uint128{uiB64, uiNonsigB0};
    }

    return return_larger_magnitude(uiA64, uiNonsigA0, uiB64, uiNonsigB0);
}

}  // namespace Intel_8086
}  // namespace internals
}  // namespace softfloat
