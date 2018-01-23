
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

#include "target.hpp"

#include "internals.hpp"
#include "softfloat/functions.h"

namespace softfloat {
namespace internals {
namespace Intel_8086 {

uint64_t
softfloat_propagateNaNF64UI(uint64_t uiA, uint64_t uiB)
{
    bool const isSigNaNA = softfloat_isSigNaNF64UI(uiA);
    bool const isSigNaNB = softfloat_isSigNaNF64UI(uiB);
    /* Make NaNs non-signaling. */
    uint64_t const uiNonsigA = uiA | UINT64_C(0x0008000000000000);
    uint64_t const uiNonsigB = uiB | UINT64_C(0x0008000000000000);

    if (isSigNaNA | isSigNaNB) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        if (isSigNaNA) {
            if (isSigNaNB) {
                goto returnLargerMag;
            }
            return isNaNF64UI(uiB) ? uiNonsigB : uiNonsigA;
        } else {
            return isNaNF64UI(uiA) ? uiNonsigA : uiNonsigB;
        }
    }
returnLargerMag:
    {
        uint64_t const uiMagA = uiNonsigA & UINT64_C(0x7FFFFFFFFFFFFFFF);
        uint64_t const uiMagB = uiNonsigB & UINT64_C(0x7FFFFFFFFFFFFFFF);
        return
            uiMagA < uiMagB ? uiNonsigB :
            uiMagB < uiMagA ? uiNonsigA :
            uiNonsigA < uiNonsigB ? uiNonsigA : uiNonsigB;
    }
}

}  // namespace Intel_8086
}  // namespace internals
}  // namespace softfloat
