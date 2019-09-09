
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

namespace softfloat {
namespace internals {
namespace fast_int64 {

template<>
uint64_t
round_pack_to<uint64_t>(bool const sign,
              uint64_t sig,
              uint64_t const sigExtra,
              softfloat_round_mode const roundingMode,
              bool const exact)
{
    bool const roundNearEven = softfloat_round_near_even == roundingMode;
    bool const doIncrement =
        !roundNearEven && softfloat_round_near_maxMag != roundingMode ?
        (sign ? softfloat_round_min : softfloat_round_max) == roundingMode && 0 != sigExtra :
        UINT64_C(0x8000000000000000) <= sigExtra;

    if (doIncrement) {
        ++sig;

        if (0 == sig) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return sign ? ui64_fromNegOverflow : ui64_fromPosOverflow;
        }

        sig &= ~static_cast<uint64_t>(!(sigExtra & UINT64_C(0x7FFFFFFFFFFFFFFF)) & roundNearEven);
    }

    if (sign && 0 != sig) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return sign ? ui64_fromNegOverflow : ui64_fromPosOverflow;
    }

    if (exact && sigExtra) {
        softfloat_raiseFlags(softfloat_flag_inexact);
    }

    return sig;
}

template
uint64_t
round_pack_to<uint64_t>(bool const sign,
                      uint64_t sig,
                      uint64_t sigExtra,
                      softfloat_round_mode roundingMode,
                      bool exact);
}  // namespace fast_int64
}  // namespace internals
}  // namespace softfloat
