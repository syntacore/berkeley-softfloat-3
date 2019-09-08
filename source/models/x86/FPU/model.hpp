/** @file

This C header file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

@copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All rights reserved.
*/
/*
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

#ifndef INTEL8086_NON_SSE_TARGET_HPP_
#define INTEL8086_NON_SSE_TARGET_HPP_

#include "Intel8086.hpp"

namespace softfloat {
namespace internals {
inline namespace x86 {
inline namespace FPU {

uint16_t
propagate_NaN(uint16_t const uiA,
              uint16_t const uiB);

/**
Interpreting `uiA' and `uiB' as the bit patterns of two 32-bit floating-
point values, at least one of which is a NaN, returns the bit pattern of
the combined NaN result.  If either `uiA' or `uiB' has the pattern of a
signaling NaN, the invalid exception is raised.
*/
inline uint32_t
propagate_NaN(uint32_t const& uiA,
              uint32_t const& uiB)
{
    static uint32_t const quietNaN_bit = UINT32_C(0x00400000);
    bool const isSigNaNA = is_sNaN(uiA);
    bool const isSigNaNB = is_sNaN(uiB);

    if (isSigNaNA || isSigNaNB) {
        softfloat_raiseFlags(softfloat_flag_invalid);

        if (isSigNaNA) {
            return quietNaN_bit | uiA;
        }
    }

    return quietNaN_bit | (is_NaN(uiA) ? uiA : uiB);
}

uint64_t
propagate_NaN(uint64_t const uiA,
              uint64_t const uiB);

template<typename Ty>
inline auto
propagate_NaN(Ty const& a,
              Ty const& b)->
    decltype(u_as_f(propagate_NaN(f_as_u(a), f_as_u(b))))
{
    return u_as_f(propagate_NaN(f_as_u(a), f_as_u(b)));
}

}  // namespace FPU
}  // namespace x86

namespace slow_int64 {
inline namespace x86 {
inline namespace FPU {
}  // namespace FPU
}  // namespace x86
}  // namespace slow_int64

namespace fast_int64 {
inline namespace x86 {
inline namespace FPU {
}  // namespace FPU
}  // namespace x86
}  // namespace fast_int64

}  // namespace internals
}  // namespace softfloat

#endif  /* INTEL8086_NON_SSE_TARGET_HPP_ */
