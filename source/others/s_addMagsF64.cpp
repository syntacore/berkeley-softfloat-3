
/** @file

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
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

#include "internals.hpp"

#include "specialize.hpp"

namespace softfloat {
namespace internals {

float64_t
softfloat_addMagsF64(uint64_t uiA, uint64_t uiB, bool signZ)
{
    using namespace softfloat::internals;
    int16_t const expA = expF64UI(uiA);
    uint64_t sigA = fracF64UI(uiA);
    int16_t const expB = expF64UI(uiB);
    uint64_t sigB = fracF64UI(uiB);
    int16_t const expDiff = expA - expB;

    if (!expDiff) {
        if (!expA) {
            return u_as_f_64(uiA + sigB);
        }

        if (expA == 0x7FF) {
            return u_as_f_64(sigA | sigB ? softfloat_propagateNaNF64UI(uiA, uiB) : uiA);
        }

        return softfloat_roundPackToF64(signZ, expA, (UINT64_C(0x0020000000000000) + sigA + sigB) << 9);
    }

    int16_t expZ;
    sigA <<= 9;
    sigB <<= 9;

    if (expDiff < 0) {
        if (expB == 0x7FF) {
            return u_as_f_64(sigB ? softfloat_propagateNaNF64UI(uiA, uiB) : packToF64UI(signZ, 0x7FF, 0));
        }

        expZ = expB;

        if (expA) {
            sigA += UINT64_C(0x2000000000000000);
        } else {
            sigA <<= 1;
        }

        sigA = softfloat_shiftRightJam64(sigA, static_cast<uint32_t>(-expDiff));
    } else {
        if (expA == 0x7FF) {
            return u_as_f_64(sigA ? softfloat_propagateNaNF64UI(uiA, uiB) : uiA);
        }

        expZ = expA;

        if (expB) {
            sigB += UINT64_C(0x2000000000000000);
        } else {
            sigB <<= 1;
        }

        sigB = softfloat_shiftRightJam64(sigB, static_cast<uint32_t>(expDiff));
    }

    uint64_t sigZ = UINT64_C(0x2000000000000000) + sigA + sigB;

    if (sigZ < UINT64_C(0x4000000000000000)) {
        --expZ;
        sigZ <<= 1;
    }

    return softfloat_roundPackToF64(signZ, expZ, sigZ);
}

}  // namespace internals
}  // namespace softfloat
