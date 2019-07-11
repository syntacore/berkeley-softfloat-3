
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

#include "target.hpp"

namespace softfloat {
namespace internals {

template<>
uint32_t
roundPackTo<uint32_t>(bool const sign,
                      uint64_t sig,
                      softfloat_round_mode const roundingMode,
                      bool const exact)
{
    uint16_t const roundIncrement =
        softfloat_round_near_even == roundingMode || softfloat_round_near_maxMag == roundingMode ? 0x800u :
        (sign ? softfloat_round_min : softfloat_round_max) == roundingMode ? 0xFFFu :
        0u;
    uint16_t const roundBits = sig & 0xFFF;
    sig += roundIncrement;

    if (0 != (UINT64_C(0xFFFFF00000000000) & sig)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return sign ? ui32_fromNegOverflow : ui32_fromPosOverflow;
    }

    bool const low_bit =
        0 == (roundBits ^ UINT16_C(0x800)) &&
        softfloat_round_near_even == roundingMode;

    uint32_t const z = (sig >> 12) & ~static_cast<uint32_t>(!!low_bit);

    if (sign && 0 != z) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return ui32_fromNegOverflow;
    }

    if (exact && 0 != roundBits) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return z;
}

template
uint32_t
roundPackTo<uint32_t>(bool const sign,
                      uint64_t sig,
                      softfloat_round_mode const roundingMode,
                      bool const exact);

}  // namespace internals
}  // namespace softfloat
