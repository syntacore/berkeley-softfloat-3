
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

#include "internals.h"

#include "specialize.h"

float32_t
softfloat_addMagsF32(uint32_t uiA, uint32_t uiB)
{
    if (isNaNF32UI(uiA) || isNaNF32UI(uiB)) {
        return u_as_f_32(softfloat_propagateNaNF32UI(uiA, uiB));
    } if (isInf32UI(uiA) || isInf32UI(uiB)) {
        return signed_inf_F32(signF32UI(uiA));
    } else {
        int16_t const expA = expF32UI(uiA);
        uint32_t const sigA = fracF32UI(uiA);
        int16_t const expB = expF32UI(uiB);
        uint32_t const sigB = fracF32UI(uiB);
        int16_t const expDiff = expA - expB;
        if (0 == expDiff) {
            /* same exponent, aligned significand fields */
            if (0 == expA) {
                /* a and b are subnormal or zero, result is simple sum, carry to exponent is valid */
                return u_as_f_32(uiA + sigB);
            } else {
                /* case aligned normalized */
                bool const signA = signF32UI(uiA);
                uint32_t const sigZ = (sigA + (UINT32_C(1) << 23)) + (sigB + (UINT32_C(1) << 23));
                if (0 == (sigZ & 1) && expA < 0xFE) {
                    /* lowest bit is 0 and sum exponent in normal range */
                    return u_as_f_32(packToF32UI(signA, expA + 1, (sigZ >> 1) & ~(~UINT32_C(0) << 23)));
                } else {
                    return softfloat_roundPackToF32(signA, expA, sigZ << 6);
                }
            }
        } else {
            /* unaligned significand fields */
            bool const signA = signF32UI(uiA);
            int16_t expZ;
            /* fractional bits before shift are [22..0] and after shift are [28..6] */
            uint32_t sigA_scaled = sigA << 6;
            uint32_t sigB_scaled = sigB << 6;
            /* bit before point is [30] */
            static uint32_t const hidden_bit = UINT32_C(1) << 29;
            if (expDiff < 0) {
                /* magnitude b greater than magnitude a */
                expZ = expB;
                /* add hidden bit and shift */
                /** @todo check suspicious  */
                sigA_scaled = softfloat_shiftRightJam32(sigA_scaled + (expA ? hidden_bit : sigA_scaled), -expDiff);
                sigB_scaled += hidden_bit;
            } else {
                /* magnitude a greater than magnitude b */
                expZ = expA;
                /* add hidden bit and shift */
                sigB_scaled = softfloat_shiftRightJam32(sigB_scaled + (expB ? hidden_bit : sigB_scaled), expDiff);
                sigA_scaled += hidden_bit;
            }
            {
                uint32_t sigZ = sigA_scaled + sigB_scaled;
                if (sigZ < 2 * hidden_bit) {
                    --expZ;
                    sigZ <<= 1;
                }
                return softfloat_roundPackToF32(signA, expZ, sigZ);
            }
        }
    }
}

