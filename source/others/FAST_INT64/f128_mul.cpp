
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

#ifndef SOFTFLOAT_FAST_INT64
#error Fast int64_t operations only
#endif

float128_t
f128_mul(float128_t const a,
         float128_t const b)
{
    using namespace softfloat::internals;
    uint64_t magBits;
    exp32_sig128 normExpSig;
    int32_t expZ;
    uint64_t sig256Z[4];
    uint64_t sigZExtra;
    uint128 sigZ;
    uint128_extra sig128Extra;
    uint128 uiZ;

    uint64_t const uiA64 = uint128(a).v64;
    uint64_t const uiA0 = uint128(a).v0;
    bool const signA = is_sign(uiA64);
    int32_t expA = expF128UI64(uiA64);

    uint128 sigA;
    sigA.v64 = fracF128UI64(uiA64);
    sigA.v0 = uiA0;

    uint64_t const uiB64 = uint128(b).v64;
    uint64_t const uiB0 = uint128(b).v0;
    bool const signB = is_sign(uiB64);
    int32_t expB = expF128UI64(uiB64);

    uint128 sigB;
    sigB.v64 = fracF128UI64(uiB64);
    sigB.v0 = uiB0;

    bool signZ = signA != signB;

    if (expA == 0x7FFF) {
        if (
            (sigA.v64 | sigA.v0) || ((expB == 0x7FFF) && (sigB.v64 | sigB.v0))
        ) {
            uiZ = softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0);
            return float128_t(uiZ);
        }

        magBits = expB | sigB.v64 | sigB.v0;
        if (!magBits) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            uiZ.v64 = defaultNaNF128UI64;
            uiZ.v0 = defaultNaNF128UI0;
            return float128_t(uiZ);
        }

        uiZ.v64 = packToF128UI64(signZ, 0x7FFF, 0);
        uiZ.v0 = 0;
        return float128_t(uiZ);
    }

    if (expB == 0x7FFF) {
        if (sigB.v64 | sigB.v0) {
            uiZ = softfloat_propagateNaNF128UI(uiA64, uiA0, uiB64, uiB0);
            return float128_t(uiZ);
        }

        magBits = expA | sigA.v64 | sigA.v0;
        if (!magBits) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            uiZ.v64 = defaultNaNF128UI64;
            uiZ.v0 = defaultNaNF128UI0;
            return float128_t(uiZ);
        }

        uiZ.v64 = packToF128UI64(signZ, 0x7FFF, 0);
        uiZ.v0 = 0;
        return float128_t(uiZ);
    }

    if (!expA) {
        if (!(sigA.v64 | sigA.v0)) {
            uiZ.v64 = packToF128UI64(signZ, 0, 0);
            uiZ.v0 = 0;
            return float128_t(uiZ);
        }

        normExpSig = softfloat_normSubnormalF128Sig(sigA.v64, sigA.v0);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (!expB) {
        if (!(sigB.v64 | sigB.v0)) {
            uiZ.v64 = packToF128UI64(signZ, 0, 0);
            uiZ.v0 = 0;
            return float128_t(uiZ);
        }

        normExpSig = softfloat_normSubnormalF128Sig(sigB.v64, sigB.v0);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    expZ = expA + expB - 0x4000;
    sigA.v64 |= UINT64_C(0x0001000000000000);
    sigB = softfloat_shortShiftLeft128(sigB, 16);
    softfloat_mul128To256M(sigA.v64, sigA.v0, sigB.v64, sigB.v0, sig256Z);
    sigZExtra = sig256Z[indexWord(4, 1)] | (sig256Z[indexWord(4, 0)] != 0);
    sigZ =
        softfloat_add128(
            sig256Z[indexWord(4, 3)], sig256Z[indexWord(4, 2)],
            sigA.v64, sigA.v0
        );

    if (UINT64_C(0x0002000000000000) <= sigZ.v64) {
        ++expZ;
        sig128Extra =
            softfloat_shortShiftRightJam128Extra(
                sigZ.v64, sigZ.v0, sigZExtra, 1);
        sigZ = sig128Extra.v;
        sigZExtra = sig128Extra.extra;
    }

    return softfloat_roundPackToF128(signZ, expZ, sigZ.v64, sigZ.v0, sigZExtra);
}

