
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

#ifdef SOFTFLOAT_FAST_INT64
#error For non-fast int64_t only
#endif

/**
@todo split to different implementations
*/
float64_t
f64_mul(float64_t a,
        float64_t b)
{
    using namespace softfloat::internals;
    bool const signA = is_sign(a);
    int16_t expA = get_exp(a);
    uint64_t sigA = get_frac(a);
    bool const signB = is_sign(b);
    int16_t expB = get_exp(b);
    uint64_t sigB = get_frac(b);
    bool const signZ = signA != signB;

    if (is_NaN(a) || is_NaN(b)) {
        return u_as_f(propagate_NaN(f_as_u(a), f_as_u(b)));
    }

    if ((is_inf(a) && is_zero(b)) || (is_zero(a) && is_inf(b))) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f(defaultNaNF64UI);
    }

    if (is_inf(a) || is_inf(b)) {
        return u_as_f(packToF64UI(signZ, 0x7FF, 0));
    }

    if (is_zero(a) || is_zero(b)) {
        return u_as_f(packToF64UI(signZ, 0, 0));
    }

    if (0 == expA) {
        exp16_sig64 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    if (0 == expB) {
        exp16_sig64 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    int16_t expZ = expA + expB - 0x3FF;
    sigA = (sigA | UINT64_C(0x0010000000000000)) << 10;
    sigB = (sigB | UINT64_C(0x0010000000000000)) << 11;
    uint32_t sig128Z[4];
    softfloat_mul64To128M(sigA, sigB, sig128Z);
    uint64_t sigZ =
        static_cast<uint64_t>(sig128Z[indexWord(4, 3)]) << 32 | sig128Z[indexWord(4, 2)];

    if (sig128Z[indexWord(4, 1)] || sig128Z[indexWord(4, 0)]) {
        sigZ |= 1;
    }

    if (sigZ < UINT64_C(0x4000000000000000)) {
        --expZ;
        sigZ <<= 1;
    }

    return softfloat_roundPackToF64(signZ, expZ, sigZ);
}
