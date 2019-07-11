
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

#include "model.hpp"

namespace softfloat {
namespace internals {
namespace Intel_8086 {

/**
Assuming at least one of the two 80-bit extended floating-point values
pointed to by `aSPtr' and `bSPtr' is a NaN, stores the combined NaN result
at the location pointed to by `zSPtr'.  If either original floating-point
value is a signaling NaN, the invalid exception is raised.
*/
/** @bug use extFloat80_t */
void
softfloat_propagateNaNExtF80M(extFloat80M const* const aSPtr,
                              extFloat80M const* const bSPtr,
                              extFloat80M* const zSPtr)
{
    if (!bSPtr) {
        if (extF80M_isSignalingNaN(aSPtr)) {
            softfloat_raiseFlags(softfloat_flag_invalid);
        }

        zSPtr->signExp = aSPtr->signExp;
        zSPtr->signif = aSPtr->signif | UINT64_C(0xC000000000000000);
    } else {
        bool const isSigNaNA = extF80M_isSignalingNaN(aSPtr);
        bool const isSigNaNB = extF80M_isSignalingNaN(bSPtr);

        if (isSigNaNA || isSigNaNB) {
            softfloat_raiseFlags(softfloat_flag_invalid);

            auto const p_src =
                !isSigNaNA ? (isNaNExtF80UI(aSPtr->signExp, aSPtr->signif) ? aSPtr : bSPtr) :
                !isSigNaNB ? (isNaNExtF80UI(bSPtr->signExp, bSPtr->signif) ? bSPtr : aSPtr) :
                (aSPtr->signExp & 0x7FFFu) < (bSPtr->signExp & 0x7FFFu) ? bSPtr :
                (bSPtr->signExp & 0x7FFFu) < (aSPtr->signExp & 0x7FFFu) ? aSPtr :
                aSPtr->signif < bSPtr->signif ? bSPtr :
                bSPtr->signif < aSPtr->signif ? aSPtr :
                aSPtr->signExp < bSPtr->signExp ? aSPtr : bSPtr;
            zSPtr->signExp = p_src->signExp;
            zSPtr->signif = p_src->signif | UINT64_C(0xC000000000000000);
        } else {
            auto const p_src =
                (aSPtr->signExp & 0x7FFFu) < (bSPtr->signExp & 0x7FFFu) ? bSPtr :
                (bSPtr->signExp & 0x7FFFu) < (aSPtr->signExp & 0x7FFFu) ? aSPtr :
                aSPtr->signif < bSPtr->signif ? bSPtr :
                aSPtr->signif > bSPtr->signif ? aSPtr :
                aSPtr->signExp < bSPtr->signExp ? aSPtr :
                bSPtr;
            zSPtr->signExp = p_src->signExp;
            zSPtr->signif = p_src->signif | UINT64_C(0xC000000000000000);
        }
    }
}

}  // namespace Intel_8086
}  // namespace internals
}  // namespace softfloat
