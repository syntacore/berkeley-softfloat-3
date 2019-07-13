
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
    uint128 const aa{a};
    uint128 const bb{b};

    bool const signA = is_sign(aa.v64);
    int32_t expA = expF128UI64(aa.v64);
    uint128 sigA{fracF128UI64(aa.v64), aa.v0};
    bool const signB = is_sign(bb.v64);
    int32_t expB = expF128UI64(bb.v64);
    uint128 sigB{fracF128UI64(bb.v64), bb.v0};
    bool const signZ = signA != signB;

    if (0x7FFF == expA) {
        if (
            0 != (sigA.v64 | sigA.v0) ||
            (0x7FFF == expB && 0 != (sigB.v64 | sigB.v0))
        ) {
            return float128_t(softfloat_propagateNaNF128UI(aa.v64, aa.v0, bb.v64, bb.v0));
        }

        uint64_t const magBits = expB | sigB.v64 | sigB.v0;

        if (0 == magBits) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return float128_t(uint128{defaultNaNF128UI64, defaultNaNF128UI0});
        }

        return float128_t(uint128{packToF128UI64(signZ, 0x7FFF, 0), 0});
    }

    if (0x7FFF == expB) {
        if (0 == (sigB.v64 | sigB.v0)) {
            return float128_t(softfloat_propagateNaNF128UI(aa, bb));
        }

        uint64_t const magBits = expA | sigA.v64 | sigA.v0;

        if (0 == magBits) {
            softfloat_raiseFlags(softfloat_flag_invalid);
            return float128_t(uint128{defaultNaNF128UI64, defaultNaNF128UI0});
        }

        return float128_t(uint128{packToF128UI64(signZ, 0x7FFF, 0), 0});
    }

    if (0 == expA) {
        if (0 == (sigA.v64 | sigA.v0)) {
            return float128_t(uint128{packToF128UI64(signZ, 0, 0), 0});
        }

        exp32_sig128 const normExpSig = softfloat_normSubnormalF128Sig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (0 == expB) {
        if (0 == (sigB.v64 | sigB.v0)) {
            return float128_t(uint128{packToF128UI64(signZ, 0, 0), 0});
        }

        exp32_sig128 const normExpSig = softfloat_normSubnormalF128Sig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int32_t expZ = expA + expB - 0x4000;
    sigA.v64 |= UINT64_C(0x0001000000000000);
    sigB = softfloat_shortShiftLeft128(sigB, 16);
    uint64_t sig256Z[4];
    softfloat_mul128To256M(sigA, sigB, sig256Z);
    uint64_t sigZExtra = sig256Z[indexWord(4, 1)] | (sig256Z[indexWord(4, 0)] != 0);
    uint128 sigZ = softfloat_add128(uint128{sig256Z[indexWord(4, 3)], sig256Z[indexWord(4, 2)]}, sigA);

    if (UINT64_C(0x0002000000000000) <= sigZ.v64) {
        ++expZ;
        uint128_extra const sig128Extra = softfloat_shortShiftRightJam128Extra(sigZ, sigZExtra, 1);
        sigZ = sig128Extra.v;
        sigZExtra = sig128Extra.extra;
    }

    return softfloat_roundPackToF128(signZ, expZ, sigZ.v64, sigZ.v0, sigZExtra);
}
