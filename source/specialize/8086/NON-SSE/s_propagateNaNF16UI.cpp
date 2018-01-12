
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

#include <stdbool.h>
#include <stdint.h>

#include "internals.hpp"
#include "specialize.hpp"
#include "softfloat/functions.h"

/**
Interpreting `uiA' and `uiB' as the bit patterns of two 16-bit floating-
point values, at least one of which is a NaN, returns the bit pattern of
the combined NaN result.  If either `uiA' or `uiB' has the pattern of a
signaling NaN, the invalid exception is raised.
*/
uint16_t
softfloat_propagateNaNF16UI(uint16_t uiA, uint16_t uiB)
{
    bool const isSigNaNA = softfloat_isSigNaNF16UI(uiA);
    bool const isSigNaNB = softfloat_isSigNaNF16UI(uiB);
    /* Make NaNs non-signaling.*/
    uint16_t const uiNonsigA = uiA | 0x0200;
    uint16_t const uiNonsigB = uiB | 0x0200;

    if (isSigNaNA | isSigNaNB) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        if (isSigNaNA) {
            if (!isSigNaNB) {
                return isNaNF16UI(uiB) ? uiNonsigB : uiNonsigA;
            }
        } else {
            return isNaNF16UI(uiA) ? uiNonsigA : uiNonsigB;
        }
    }
    {
        uint16_t const uiMagA = uiNonsigA & 0x7FFF;
        uint16_t const uiMagB = uiNonsigB & 0x7FFF;
        return
            uiMagA < uiMagB ? uiNonsigB :
            uiMagB < uiMagA ? uiNonsigA :
            uiNonsigA < uiNonsigB ? uiNonsigA : uiNonsigB;
    }
}
