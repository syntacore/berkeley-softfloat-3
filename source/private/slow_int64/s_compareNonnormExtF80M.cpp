
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
namespace slow_int64 {

/**
@bug use extFloat80_t
*/
int
compare_non_norm_M_extF80(extFloat80M const* const aSPtr,
                          extFloat80M const* const bSPtr)
{
    uint64_t sigA = aSPtr->signif;
    bool signB = is_sign(bSPtr->signExp);
    uint64_t sigB = bSPtr->signif;

    if (0 != ((aSPtr->signExp ^ bSPtr->signExp) & 0x8000)) {
        if (0 == (sigA | sigB)) {
            return 0;
        }

        return signB ? 1 : -1;
    }

    int32_t expA = get_exp(*aSPtr);
    int32_t expB = get_exp(*bSPtr);

    if (0x7FFF == expA) {
        if (0x7FFF == expB) {
            return 0;
        }

        return signB ? -1 : 1;
    }

    if (0x7FFF == expB) {
        return signB ? 1 : -1;
    }

    if (0 == expA) {
        expA = 1;
    }

    if (0 == (sigA & UINT64_C(0x8000000000000000))) {
        expA = 0 != sigA ? expA + norm_M_extF80Sig(&sigA) : -128;
    }

    if (0 == expB) {
        expB = 1;
    }

    if (0 == (sigB & UINT64_C(0x8000000000000000))) {
        expB = 0 != sigB ? expB + norm_M_extF80Sig(&sigB) : -128;
    }

    if (signB) {
        if (expA < expB) {
            return 1;
        }

        if (expB < expA || sigB < sigA) {
            return -1;
        }
    } else {
        if (expB < expA) {
            return 1;
        }

        if (expA < expB || sigA < sigB) {
            return -1;
        }
    }

    return !!(sigA != sigB);
}

}  // namespace slow_int64
}  // namespace internals
}  // namespace softfloat
