
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

#if (SOFTFLOAT_FAST_INT64)

float64_t
f64_rem(float64_t const a,
        float64_t const b)
{
    using namespace softfloat::internals::fast_int64;

    bool const signA = is_sign(a);
    int16_t expA = get_exp(a);
    uint64_t sigA = get_frac(a);
    int16_t expB = get_exp(b);
    uint64_t sigB = get_frac(b);

    if (is_NaN(a) || is_NaN(b)) {
        return propagate_NaN(a, b);
    }

    if (is_inf(a)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f(defaultNaNF64UI);
    }

    if (is_inf(b)) {
        return a;
    }

    if (expA < expB - 1) {
        return a;
    }

    if (is_zero(b)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f(defaultNaNF64UI);
    }

    if (is_zero(a)) {
        return a;
    }

    if (0 == expB) {
        exp16_sig64 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    if (0 == expA) {
        exp16_sig64 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    uint64_t rem = sigA | UINT64_C(0x0010000000000000);
    sigB |= UINT64_C(0x0010000000000000);
    int16_t expDiff = expA - expB;
    uint32_t q;
    uint64_t altRem;

    if (expDiff < 1) {
        if (expDiff < -1) {
            return a;
        }

        sigB <<= 9;

        if (0 != expDiff) {
            rem <<= 8;
            q = 0;
        } else {
            rem <<= 9;
            q = 0u + !!(sigB <= rem);

            if (q) {
                rem -= sigB;
            }
        }
    } else {
        uint32_t const recip32 = approxRecip32_1(static_cast<uint32_t>(sigB >> 21));
        /*
        Changing the shift of `rem' here requires also changing the initial
        subtraction from `expDiff'.
        */
        rem <<= 9;
        expDiff -= 30;
        /*
        The scale of `sigB' affects how many bits are obtained during each
        cycle of the loop.  Currently this is 29 bits per loop iteration,
        the maximum possible.
        */
        sigB <<= 9;

        uint64_t q64;

        for (;;) {
            q64 = static_cast<uint32_t>(rem >> 32) * static_cast<uint64_t>(recip32);

            if (expDiff < 0) {
                break;
            }

            q = (q64 + 0x80000000) >> 32;
            rem <<= 29;
            rem -= q * static_cast<uint64_t>(sigB);

            if (rem & UINT64_C(0x8000000000000000)) {
                rem += sigB;
            }

            expDiff -= 29;
        }

        /* `expDiff' cannot be less than -29 here. */
        q = static_cast<uint32_t>(q64 >> 32) >> (~expDiff & 31);
        rem = (rem << (expDiff + 30)) - q * static_cast<uint64_t>(sigB);

        if (0 != (rem & UINT64_C(0x8000000000000000))) {
            altRem = rem + sigB;
            uint64_t const meanRem = rem + altRem;

            if ((meanRem & UINT64_C(0x8000000000000000)) || (!meanRem && (q & 1))) {
                rem = altRem;
            }

            bool signRem = signA;

            if (0 != (rem & UINT64_C(0x8000000000000000))) {
                signRem = !signRem;
                rem = static_cast<uint64_t>(-static_cast<int64_t>(rem));
            }

            return normRoundPackToF64(signRem, expB, rem);
        }
    }

    do {
        altRem = rem;
        ++q;
        rem -= sigB;
    } while (0 == (rem & UINT64_C(0x8000000000000000)));


    uint64_t const meanRem = rem + altRem;

    if (0 != (meanRem & UINT64_C(0x8000000000000000)) || (0 == meanRem && 0 != (q & 1))) {
        rem = altRem;
    }

    bool signRem = signA;

    if (rem & UINT64_C(0x8000000000000000)) {
        signRem = !signRem;
        rem = static_cast<uint64_t>(-static_cast<int64_t>(rem));
    }

    return normRoundPackToF64(signRem, expB, rem);
}

#else

using namespace softfloat::internals;

namespace {
static float64_t
selectRem(bool const signA,
          int16_t expB,
          uint64_t rem,
          uint64_t altRem,
          uint32_t q)
{
    uint64_t const meanRem = rem + altRem;

    if (0 != (meanRem & UINT64_C(0x8000000000000000)) || (0 == meanRem && 0 != (q & 1))) {
        rem = altRem;
    }

    bool signRem = signA;

    if (0 != (rem & UINT64_C(0x8000000000000000))) {
        signRem = !signRem;
        rem = static_cast<uint64_t>(-static_cast<int64_t>(rem));
    }

    return normRoundPackToF64(signRem, expB, rem);
}

}  // namespace

float64_t
f64_rem(float64_t const a,
        float64_t const b)
{
    if (is_NaN(a) || is_NaN(b)) {
        return propagate_NaN(a, b);
    }

    if (is_inf(a)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f(defaultNaNF64UI);
    }

    if (is_inf(b)) {
        return a;
    }

    if (is_zero(b)) {
        softfloat_raiseFlags(softfloat_flag_invalid);
        return u_as_f(defaultNaNF64UI);
    }

    bool const signA = is_sign(a);
    int16_t expA = get_exp(a);
    uint64_t sigA = get_frac(a);
    int16_t expB = get_exp(b);
    uint64_t sigB = get_frac(b);

    if (expA < expB - 1) {
        return a;
    }

    if (is_zero(a)) {
        return a;
    }

    if (0 == expB) {
        exp16_sig64 const normExpSig(sigB);
        expB = normExpSig.exp;
        sigB = normExpSig.sig;
    }

    if (0 == expA) {
        exp16_sig64 const normExpSig(sigA);
        expA = normExpSig.exp;
        sigA = normExpSig.sig;
    }

    uint64_t rem = sigA | UINT64_C(0x0010000000000000);
    sigB |= UINT64_C(0x0010000000000000);
    int16_t expDiff = expA - expB;
    uint32_t q;

    if (expDiff < 1) {
        if (expDiff < -1) {
            return a;
        }

        sigB <<= 9;

        if (expDiff) {
            rem <<= 8;
            q = 0;
        } else {
            rem <<= 9;
            q = 0u + !!(sigB <= rem);

            if (q) {
                rem -= sigB;
            }
        }
    } else {
        uint32_t const recip32 = approxRecip32_1(static_cast<uint32_t>(sigB >> 21));
        /*
        Changing the shift of `rem' here requires also changing the initial
        subtraction from `expDiff'.
        */
        rem <<= 9;
        expDiff -= 30;
        /*
        The scale of `sigB' affects how many bits are obtained during each
        cycle of the loop.  Currently this is 29 bits per loop iteration,
        the maximum possible.
        */
        sigB <<= 9;
        uint64_t q64;

        for (;;) {
            q64 = static_cast<uint32_t>(rem >> 32) * static_cast<uint64_t>(recip32);

            if (expDiff < 0) {
                break;
            }

            q = (q64 + 0x80000000) >> 32;
            rem = static_cast<uint64_t>(static_cast<uint32_t>(rem >> 3)) << 32;
            rem -= q * static_cast<uint64_t>(sigB);

            if (0 != (rem & UINT64_C(0x8000000000000000))) {
                rem += sigB;
            }

            expDiff -= 29;
        }

        /* `expDiff' cannot be less than -29 here. */
        q = static_cast<uint32_t>(q64 >> 32) >> (~expDiff & 31);
        rem = (rem << (expDiff + 30)) - q * static_cast<uint64_t>(sigB);

        if (0 != (rem & UINT64_C(0x8000000000000000))) {
            return selectRem(signA, expB, rem, rem + sigB, q);
        }
    }

    uint64_t altRem;

    do {
        altRem = rem;
        ++q;
        rem -= sigB;
    } while (0 == (rem & UINT64_C(0x8000000000000000)));

    return selectRem(signA, expB, rem, altRem, q);

}

#endif
