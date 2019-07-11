
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

#include "model.hpp"

/** @todo split implementations */
uint64_t
f64_to_ui64(float64_t const a,
            softfloat_round_mode const roundingMode,
            bool const exact)
{
    using namespace softfloat::internals;
    uint64_t const uiA = f_as_u_64(a);
    bool const sign = is_sign(uiA);
    int16_t const exp = expF64UI(uiA);
    uint64_t const sig = fracF64UI(uiA) | (0 != exp ? UINT64_C(0x0010000000000000) : UINT64_C(0));
    int16_t const shiftDist = 0x433 - exp;

    if (shiftDist <= 0) {
        if (shiftDist < -11) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return
                exp == 0x7FF && fracF64UI(uiA) ? ui64_fromNaN :
                sign ? ui64_fromNegOverflow : ui64_fromPosOverflow;
        }

        auto const sig_1 = sig << -shiftDist;
        uint32_t extSig[3];
        extSig[indexWord(3, 0)] = 0;
        extSig[indexWord(3, 2)] = sig_1 >> 32;
        extSig[indexWord(3, 1)] = static_cast<uint32_t>(sig_1);
        return roundPackMTo<uint64_t>(sign, extSig, roundingMode, exact);
    }

    uint32_t extSig[3];
    extSig[indexWord(3, 0)] = 0;
    extSig[indexWord(3, 2)] = sig >> 32;
    extSig[indexWord(3, 1)] = static_cast<uint32_t>(sig);
    softfloat_shiftRightJam96M(extSig, static_cast<uint8_t>(shiftDist), extSig);
    return roundPackMTo<uint64_t>(sign, extSig, roundingMode, exact);
}
