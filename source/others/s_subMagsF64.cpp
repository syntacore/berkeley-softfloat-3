
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
#include "softfloat/functions.h"

float64_t
softfloat_subMagsF64(uint64_t uiA, uint64_t uiB, bool signZ)
{
    int16_t expA = expF64UI(uiA);
    uint64_t const sigA = fracF64UI(uiA);
    int16_t expB = expF64UI(uiB);
    uint64_t const sigB = fracF64UI(uiB);

    int16_t const expDiff = expA - expB;
    if (0 == expDiff) {
        if (expA == 0x7FF) {
            if (0 != sigA || 0 != sigB) {
                return u_as_f_64(softfloat_propagateNaNF64UI(uiA, uiB));
            } else {
                softfloat_raiseFlags(softfloat_flag_invalid);
                return u_as_f_64(defaultNaNF64UI);
            }
        }
        int64_t sigDiff = sigA - sigB;
        if (0 == sigDiff) {
            return u_as_f_64(packToF64UI(softfloat_round_min == softfloat_roundingMode, 0, 0));
        } else {
            if (0 != expA) {
                --expA;
            }
            if (sigDiff < 0) {
                signZ = !signZ;
                sigDiff = -sigDiff;
            }
            int8_t const shiftDist = softfloat_countLeadingZeros64(sigDiff) - 11;
            int16_t const expZ = expA - shiftDist;
            return
                u_as_f_64(expZ < 0 ?
                          packToF64UI(signZ, 0, sigDiff << expA) :
                          packToF64UI(signZ, expZ, sigDiff << shiftDist));
        }
    } else {
        uint64_t sigA_shifted = sigA << 10;
        uint64_t sigB_shifted = sigB << 10;
        if (expDiff < 0) {
            signZ = !signZ;
            if (0x7FF == expB) {
                return
                    u_as_f_64(0 != sigB_shifted ?
                              softfloat_propagateNaNF64UI(uiA, uiB) :
                              packToF64UI(signZ, 0x7FF, 0));
            } else {
                return
                    softfloat_normRoundPackToF64(signZ,
                                                 expB - 1,
                                                 (sigB_shifted | UINT64_C(0x4000000000000000)) -
                                                 softfloat_shiftRightJam64(sigA_shifted + (expA ? UINT64_C(0x4000000000000000) : sigA_shifted), -expDiff));
            }
        } else {
            if (0x7FF == expA) {
                return u_as_f_64(0 == sigA_shifted ? uiA : softfloat_propagateNaNF64UI(uiA, uiB));
            } else {
                return
                    softfloat_normRoundPackToF64(signZ,
                                                 expA - 1,
                                                 (sigA_shifted | UINT64_C(0x4000000000000000)) -
                                                 softfloat_shiftRightJam64(sigB_shifted + (expB ? UINT64_C(0x4000000000000000) : sigB_shifted), expDiff));
            }
        }
    }
}
