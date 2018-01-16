
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

#include "primitives/functions.hpp"
#include "specialize.hpp"
#include "softfloat/functions.h"

namespace softfloat {
namespace Intel_8086 {

/**
Assuming the unsigned integer formed from concatenating `uiA64' and `uiA0'
has the bit pattern of a 128-bit floating-point NaN, converts this NaN to
the common NaN form, and stores the resulting common NaN at the location
pointed to by `zPtr'.  If the NaN is a signaling NaN, the invalid exception
is raised.
*/
commonNaN
softfloat_f128UIToCommonNaN(uint64_t uiA64,
                            uint64_t uiA0)
{
    if (softfloat_isSigNaNF128UI(uiA64, uiA0)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }
    uint128 const NaNSig = softfloat_shortShiftLeft128(uiA64, uiA0, 16);
    commonNaN z;
    z.sign = 0 != (uiA64 >> 63);
    z.v64 = NaNSig.v64;
    z.v0 = NaNSig.v0;
    return z;
}

}  // namespace Intel_8086
}  // namespace softfloat
