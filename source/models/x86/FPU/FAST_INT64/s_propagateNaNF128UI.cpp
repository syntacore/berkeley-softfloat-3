
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

namespace softfloat {
namespace internals {
namespace Intel_8086 {

namespace {
static uint64_t const non_signal_bit = UINT64_C(0x0000'8000'0000'0000);
static uint64_t const suppress_sign_mask = UINT64_C(0x7FFF'FFFF'FFFF'FFFF);
}

uint128
softfloat_propagateNaNF128UI(uint64_t const& uiA64,
                             uint64_t const& uiA0,
                             uint64_t const& uiB64,
                             uint64_t const& uiB0)
{
    bool const isSigNaNA = softfloat_isSigNaNF128UI(uiA64, uiA0);
    bool const isSigNaNB = softfloat_isSigNaNF128UI(uiB64, uiB0);
    /* Make NaNs non-signaling. */
    uint64_t const uiNonsigA64 = uiA64 | non_signal_bit;
    uint64_t const uiNonsigB64 = uiB64 | non_signal_bit;
    uint128 nsA(uiNonsigA64, uiA0);
    uint128 nsB(uiNonsigB64, uiB0);

    if (isSigNaNA || isSigNaNB) {
        softfloat_raiseFlags(softfloat_flag_invalid);

        if (!isSigNaNA) {
            return isNaNF128UI(uiA64, uiA0) ? nsA : nsB;
        }

        assert(isSigNaNA);

        if (!isSigNaNB) {
            return isNaNF128UI(uiB64, uiB0) ? nsB : nsA;
        }

        assert(isSigNaNA && isSigNaNB);
    }

    uint64_t const uiMagA64 = uiNonsigA64 & suppress_sign_mask;
    uint64_t const uiMagB64 = uiNonsigB64 & suppress_sign_mask;

    if (uiMagA64 < uiMagB64) {
        return nsB;
    }

    assert(uiMagB64 <= uiMagA64);

    if (uiMagB64 < uiMagA64) {
        return nsA;
    }

    assert(uiMagB64 == uiMagA64);

    if (uiA0 < uiB0) {
        return nsB;
    }

    assert(uiB0 <= uiA0);

    if (uiB0 < uiA0) {
        return nsA;
    }

    assert(uiB0 == uiA0);

    return
        uiNonsigA64 < uiNonsigB64 ? nsA : nsB;
}

}  // namespace Intel_8086
}  // namespace internals
}  // namespace softfloat
