
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015 The Regents of the University of
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

#include "primitives/types.hpp"
#include "specialize.hpp"
#include "softfloat/functions.h"

namespace softfloat {
namespace internals {
namespace riscv {

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
                              extFloat80M* zSPtr)
{
    uint16_t ui64 = aSPtr->signExp;
    uint64_t ui0 = aSPtr->signif;

    if (softfloat_isSigNaNExtF80UI(ui64, ui0) || (bSPtr && (ui64 = bSPtr->signExp, ui0 = bSPtr->signif, softfloat_isSigNaNExtF80UI(ui64, ui0)))) {
        softfloat_raiseFlags(softfloat_flag_invalid);
    }

    zSPtr->signExp = defaultNaNExtF80UI64;
    zSPtr->signif = defaultNaNExtF80UI0;
}

}  // namespace riscv 
}  // namespace internals
}  // namespace softfloat
