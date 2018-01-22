
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

#include "specialize.hpp"

#include "internals.hpp"
#include "softfloat/functions.h"

namespace softfloat {
namespace internals {
namespace Intel_8086 {

uint32_t
softfloat_propagateNaNF32UI(uint32_t uiA, uint32_t uiB)
{
    bool const isSigNaNA = softfloat_isSigNaNF32UI(uiA);
    bool const isSigNaNB = softfloat_isSigNaNF32UI(uiB);
    /* Make NaNs non-signaling. */
    uint32_t const uiNonsigA = uiA | 0x00400000;
    uint32_t const uiNonsigB = uiB | 0x00400000;

    if (isSigNaNA || isSigNaNB) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        if (!isSigNaNA) {
            return isNaNF32UI(uiA) ? uiNonsigA : uiNonsigB;
        } else if (!isSigNaNB) {
            return isNaNF32UI(uiB) ? uiNonsigB : uiNonsigA;
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

}  // namespace Intel_8086
}  // namespace internals
}  // namespace softfloat
