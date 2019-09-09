/** @file

This C header file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

@copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
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
*/
/*
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

#ifndef SOFTFLOAT_INTERNALS_HPP_
#define SOFTFLOAT_INTERNALS_HPP_

#include "softfloat/functions.h"
#include "primitives/functions.hpp"

#include <cassert>
#include <type_traits>

#if !defined(SOFTFLOAT_FAST_INT64) || !defined(SOFTFLOAT_FAST_DIV64TO32) || !defined(SOFTFLOAT_FAST_DIV32TO16)
#error SOFTFLOAT_FAST_INT64, SOFTFLOAT_FAST_DIV64TO32 and SOFTFLOAT_FAST_DIV32TO16 should be defined as 0/1
#endif

#ifndef THREAD_LOCAL
#define THREAD_LOCAL thread_local
#endif

namespace softfloat {
namespace internals {

enum Mul_add_operations {
    softfloat_mulAdd_madd = 0,
    softfloat_mulAdd_subC = 1,
    softfloat_mulAdd_subProd = 2
};

struct exp8_sig16 {
    exp8_sig16() = default;

    constexpr
    exp8_sig16(int8_t a_exp,
               uint16_t a_sig)
        : exp(a_exp)
        , sig(a_sig)
    {}

    explicit
    exp8_sig16(uint16_t const sig)
        : exp(static_cast<int8_t>(1 - (count_leading_zeros(sig) - 5)))
        , sig(static_cast<uint16_t>(sig << (count_leading_zeros(sig) - 5)))
    {}

    int8_t exp;
    uint16_t sig;
};

struct exp16_sig32 {
    exp16_sig32() = default;

    constexpr
    exp16_sig32(int16_t a_exp,
                uint32_t a_sig)
        : exp(a_exp)
        , sig(a_sig)
    {}

    explicit
    exp16_sig32(uint32_t const sig)
        : exp{static_cast<int16_t>(1 - (count_leading_zeros(sig) - 8))}
        , sig{sig << (count_leading_zeros(sig) - 8)}
    {}

    int16_t exp;
    uint32_t sig;
};

struct exp16_sig64 {
    exp16_sig64() = default;
    constexpr exp16_sig64(int16_t a_exp,
                          uint64_t a_sig)
        : exp(a_exp)
        , sig(a_sig)
    {}

    explicit exp16_sig64(uint64_t const sig)
        : exp{static_cast<int16_t>(1 - (count_leading_zeros(sig) - 11))}
        , sig{sig << (count_leading_zeros(sig) - 11)}
    {}

    int16_t exp;
    uint64_t sig;
};

/**
Rounding precision for 80-bit extended double-precision floating-point.
Valid values are 32, 64, and 80.
*/
extern /*THREAD_LOCAL*/ uint8_t extF80_rounding_precision;
extern /*THREAD_LOCAL*/ softfloat_tininess detect_tininess;

template<typename Ty>
inline constexpr auto
f_as_u(Ty const v)->decltype(v.v)
{
    return v.v;
}

inline constexpr float16_t
u_as_f(uint16_t v)
{
    return float16_t{v};
}

inline constexpr float32_t
u_as_f(uint32_t v)
{
    return float32_t{v};
}

inline constexpr float64_t
u_as_f(uint64_t v)
{
    return float64_t{v};
}

template<typename Ty>
inline constexpr
typename std::enable_if<std::is_integral<Ty>::value, bool>::type
is_sign(Ty const& v)
{
    return static_cast<typename std::make_signed<Ty>::type>(v) < 0;
}

template<typename Ty>
inline constexpr auto
is_sign(Ty const v)->decltype(is_sign(f_as_u(v)))
{
    return is_sign(f_as_u(v));
}

template<typename Ty>
Ty
round_pack_to(bool,
              uint64_t,
              softfloat_round_mode,
              bool);

float16_t
round_pack_to_F16(bool,
                  int16_t,
                  uint16_t);
float16_t
norm_round_pack_to_F16(bool,
                       int16_t,
                       uint16_t const&);

float16_t
add_magnitudes(uint16_t const&,
               uint16_t const&);
float16_t
sub_magnitudes(uint16_t const&,
               uint16_t const&);

float32_t
round_pack_to_F32(bool,
                  int16_t,
                  uint32_t);

float64_t
round_pack_to_F64(bool,
                  int16_t,
                  uint64_t);

namespace fast_int64 {

struct exp32_sig64 {
    exp32_sig64() = default;

    constexpr
    exp32_sig64(int32_t a_exp,
                uint64_t a_sig)
        : exp(a_exp)
        , sig(a_sig)
    {}

    int32_t exp;
    uint64_t sig;
};

struct exp32_sig128 {
    exp32_sig128() = default;

    constexpr
    exp32_sig128(int32_t a_exp,
                 uint128 const& a_sig)
        : exp(a_exp)
        , sig(a_sig)
    {}

    int32_t exp;
    uint128 sig;
};

template<typename Ty>
Ty
round_pack_to(bool,
              uint64_t,
              uint64_t,
              softfloat_round_mode,
              bool);

template<>
int64_t
round_pack_to<int64_t>(bool const sign,
                       uint64_t sig,
                       uint64_t sigExtra,
                       softfloat_round_mode roundingMode,
                       bool exact);
template<>
uint64_t
round_pack_to<uint64_t>(bool const sign,
                        uint64_t sig,
                        uint64_t sigExtra,
                        softfloat_round_mode roundingMode,
                        bool exact);

inline
exp32_sig64
norm_subnormal_extF80Sig(uint64_t const& sig)
{
    auto const shiftDist = count_leading_zeros(sig);
    return exp32_sig64{-shiftDist, sig << shiftDist};
}

extFloat80_t
round_pack_to_extF80(bool,
                     int32_t,
                     uint64_t,
                     uint64_t,
                     uint8_t);
extFloat80_t
norm_round_pack_to_extF80(bool,
                          int32_t,
                          uint64_t,
                          uint64_t,
                          uint8_t);

extFloat80_t
add_magnitudes(uint16_t,
               uint64_t,
               uint16_t,
               uint64_t,
               bool);
extFloat80_t
sub_magnitudes(uint16_t,
               uint64_t,
               uint16_t,
               uint64_t,
               bool);

exp32_sig128
norm_subnormal_F128Sig(uint64_t,
                       uint64_t);

inline
exp32_sig128
norm_subnormal_F128Sig(uint128 const& a)
{
    return norm_subnormal_F128Sig(a.v64, a.v0);
}

float128_t
round_pack_to_F128(bool,
                   int32_t,
                   uint64_t,
                   uint64_t,
                   uint64_t);
float128_t
norm_round_pack_to_F128(bool,
                        int32_t,
                        uint64_t,
                        uint64_t);

float128_t
add_magnitudes(uint64_t,
               uint64_t,
               uint64_t,
               uint64_t,
               bool);

inline
float128_t
add_magnitudes(uint128 const& a,
               uint128 const& b,
               bool const signZ)
{
    return add_magnitudes(a.v64, a.v0, b.v64, b.v0, signZ);
}

float128_t
sub_magnitudes(uint64_t,
               uint64_t,
               uint64_t,
               uint64_t,
               bool);
inline
float128_t
sub_magnitudes(uint128 const& a,
               uint128 const& b,
               bool const signZ)
{
    return sub_magnitudes(a.v64, a.v0, b.v64, b.v0, signZ);
}

float128_t
mul_add_F128(Mul_add_operations const op,
             uint64_t const& uiA64,
             uint64_t const& uiA0,
             uint64_t const& uiB64,
             uint64_t const& uiB0,
             uint64_t const& uiC64,
             uint64_t const& uiC0);

inline constexpr int32_t
exp_F128_UI64(uint64_t const& a64)
{
    return static_cast<int32_t>(a64 >> 48) & INT32_C(0x7FFF);
}

inline constexpr uint64_t
frac_F128_UI64(uint64_t const& a64)
{
    return a64 & UINT64_C(0x0000FFFFFFFFFFFF);
}

inline constexpr uint64_t
pack_to_F128_UI64(bool const sign,
                  int32_t const& exp,
                  uint64_t const& sig64)
{
    return
        (static_cast<uint64_t>(sign) << 63) +
        (static_cast<uint64_t>(exp) << 48) +
        sig64;
}

inline constexpr bool
is_NaN(uint64_t const& a64,
       uint64_t const& a0)
{
    return
        0 == (~a64 & UINT64_C(0x7FFF000000000000)) &&
        (0 != a0 || 0 != (a64 & UINT64_C(0x0000FFFFFFFFFFFF)));
}

inline constexpr bool
is_NaN(uint128 const& a)
{
    return is_NaN(a.v64, a.v0);
}

}  // namespace fast_int64

namespace slow_int64 {

template<typename Ty>
Ty
round_pack_to_M(bool,
                uint32_t const*,
                softfloat_round_mode,
                bool);

template<>
int64_t
round_pack_to_M<int64_t>(bool const sign,
                         uint32_t const* const extSigPtr,
                         softfloat_round_mode const roundingMode,
                         bool const exact);

template<>
uint64_t
round_pack_to_M<uint64_t>(bool const sign,
                          uint32_t const* const extSigPtr,
                          softfloat_round_mode const roundingMode,
                          bool const exact);

bool
try_propagate_NaN_M_extF80(extFloat80M const*,
                           extFloat80M const*,
                           extFloat80M*);
void
invalid_M_extF80(extFloat80M*);

int
norm_M_extF80Sig(uint64_t*);

void
round_pack_to_M_extF80M(bool,
                        int32_t,
                        uint32_t* const,
                        uint8_t const&,
                        extFloat80M* const);
void
norm_round_pack_to_M_extF80(bool,
                            int32_t,
                            uint32_t*,
                            uint8_t const&,
                            extFloat80M*);

void
add_M_extF80(extFloat80M const* const,
             extFloat80M const* const,
             extFloat80M* const,
             bool);

int
compare_non_norm_M_extF80(extFloat80M const*,
                          extFloat80M const*);

bool is_NaN_M_F128(uint32_t const*);

bool
try_propagate_NaN_M_F128(uint32_t const*,
                         uint32_t const*,
                         uint32_t*);

void
invalid_M_F128(uint32_t*);

int
shift_norm_sig_M_F128(uint32_t const*,
                      uint8_t,
                      uint32_t*);

void
round_pack_to_M_F128(bool,
                     int32_t,
                     uint32_t*,
                     uint32_t*);
void
norm_round_pack_to_M_F128(bool,
                          int32_t,
                          uint32_t*,
                          uint32_t*);

void
add_M_F128(uint32_t const*,
           uint32_t const*,
           uint32_t*,
           bool);
void
mul_add_M_F128(Mul_add_operations,
               uint32_t const*,
               uint32_t const*,
               uint32_t const*,
               uint32_t*);

inline constexpr uint16_t
exp_F128_UI96(uint32_t const& a96)
{
    return static_cast<uint16_t>((a96 >> 16) & UINT16_C(0x7FFF));
}

inline constexpr uint32_t
frac_F128_UI96(uint32_t const& a96)
{
    return a96 & 0x0000FFFF;
}

inline constexpr uint32_t
pack_to_F128_UI96(bool const sign,
                  unsigned const expnt,
                  uint32_t const& sig96)
{
    return (static_cast<uint32_t>((!!sign) << 31) + (static_cast<uint32_t>(expnt) << 16)) + sig96;
}
}  // namespace slow_int64

inline constexpr int8_t
get_exp(uint16_t const& a)
{
    return static_cast<int8_t>((a >> 10) & 0x1F);
}

inline constexpr int16_t
get_exp(uint32_t const& a)
{
    return static_cast<int16_t>((a >> 23) & ~(~UINT32_C(0) << 8));
}

inline constexpr int16_t
get_exp(uint64_t const& a)
{
    return static_cast<int16_t>((a >> 52) & ~(~UINT32_C(0) << 11));
}

template<typename Ty>
inline constexpr auto
get_exp(Ty const a)->decltype(get_exp(f_as_u(a)))
{
    return get_exp(f_as_u(a));
}

inline constexpr uint16_t
get_frac(uint16_t const& a)
{
    return a & 0x03FFu;
}

inline constexpr uint32_t
get_frac(uint32_t const& a)
{
    return a & ~(~UINT32_C(0) << 23);
}

inline constexpr uint64_t
get_frac(uint64_t const& a)
{
    return a & ~(~UINT64_C(0) << 52);
}

inline constexpr uint16_t
exp_extF80_UI64(uint16_t const& a64)
{
    return static_cast<uint16_t>(0x7FFF & a64);
}

template<typename Ty>
inline constexpr auto
get_frac(Ty const a)->decltype(get_frac(f_as_u(a)))
{
    return get_frac(f_as_u(a));
}

inline constexpr bool
is_finite(uint16_t const& a)
{
    return 0x1F != get_exp(a);
}

inline constexpr bool
is_finite(uint32_t const& a)
{
    return 0xFF != get_exp(a);
}

inline constexpr bool
is_finite(uint64_t const& a)
{
    return 0x7FF != get_exp(a);
}

template<typename Ty>
inline constexpr auto
is_finite(Ty const a)->decltype(is_finite(f_as_u(a)))
{
    return is_finite(f_as_u(a));
}

template<typename Ty>
inline constexpr auto
is_denormalized(Ty const a)->decltype(get_exp(a) == 0)
{
    return 0 == get_exp(a);
}

template<typename Ty>
inline constexpr auto
is_zero(Ty const a)->decltype(is_denormalized(a)&& 0 == get_frac(a))
{
    return is_denormalized(a) && 0 == get_frac(a);
}

template<typename Ty>
inline constexpr auto
is_subnormal(Ty const a)->decltype(is_denormalized(a)&& 0 != get_frac(a))
{
    return is_denormalized(a) && 0 != get_frac(a);
}

template<typename Ty>
inline constexpr auto
is_inf(Ty const a)->decltype(!is_finite(a)&& 0 == get_frac(a))
{
    return !is_finite(a) && 0 == get_frac(a);
}

template<typename Ty>
inline constexpr auto
is_NaN(Ty const a)->decltype(!is_finite(a)&& 0 != get_frac(a))
{
    return !is_finite(a) && 0 != get_frac(a);
}

inline constexpr bool
is_NaN(uint16_t const& a64,
       uint64_t const& a0)
{
    return
        UINT16_C(0x7FFF) == (UINT16_C(0x7FFF) & a64) &&
        0 != (UINT64_C(0x7FFFFFFFFFFFFFFF) & a0);
}

inline constexpr bool
is_NaN(extFloat80_t const& a)
{
    return is_NaN(a.signExp, a.signif);
}

inline constexpr bool
is_sNaN(uint16_t const uiA)
{
    return
        UINT16_C(0x7C00) == (UINT16_C(0x7E00) & uiA) &&
        0 != (UINT16_C(0x01FF) & uiA);
}

inline constexpr bool
is_sNaN(uint32_t const uiA)
{
    return
        (~(~UINT32_C(0) << 8) << 23) == (((~(~UINT32_C(0) << 8) << 23) | (~(~UINT32_C(0) << 1) << 22)) & uiA) &&
        0 != (~(~UINT32_C(0) << 23) & uiA);
}

inline constexpr bool
is_sNaN(uint64_t const uiA)
{
    return
        UINT64_C(0x7FF0000000000000) == (uiA & UINT64_C(0x7FF8000000000000)) &&
        0 != (uiA & UINT64_C(0x0007FFFFFFFFFFFF));
}

inline constexpr bool
is_sNaN(uint16_t const uiA64,
        uint64_t const uiA0)
{
    return
        UINT16_C(0x7FFF) == (UINT16_C(0x7FFF) & uiA64) &&
        0 == (UINT64_C(0x4000000000000000) & uiA0) &&
        0 != (UINT64_C(0x3FFFFFFFFFFFFFFF) & uiA0);
}

inline constexpr bool
is_sNaN(extFloat80_t const& a)
{
    return is_sNaN(a.signExp, a.signif);
}

template<typename Ty>
inline constexpr auto
is_sNaN(Ty const a)->decltype(is_sNaN(f_as_u(a)))
{
    return is_sNaN(f_as_u(a));
}

inline constexpr uint16_t
pack_to_F16_UI(bool sign,
               int8_t const& expnt,
               uint16_t const& sig)
{
    return
        static_cast<uint16_t>(
            (static_cast<uint16_t>(!!sign) << 15) +
            (static_cast<uint16_t>(expnt) << 10) +
            sig);
}

inline constexpr uint32_t
pack_to_F32_UI(bool sign,
               int16_t const& expnt,
               uint32_t const& sgnf)
{
    return
        (static_cast<uint32_t>(!!sign) << 31) |
        ((static_cast<uint32_t>(expnt) << 23) + sgnf);
}

inline constexpr uint64_t
pack_to_F64_UI(bool sign, int16_t const& expnt, uint64_t const& sgnf)
{
    return
        (static_cast<uint64_t>(!!sign) << 63) |
        ((static_cast<uint64_t>(expnt) << 52) + sgnf);
}

inline constexpr uint16_t
pack_to_extF80_UI64(bool const sign,
                    uint16_t const expnt)
{
    return static_cast<uint16_t>((!!sign << 15) | expnt);
}

template<typename Ty>
Ty
make_signed_inf(bool sign);

template<>
inline constexpr float16_t
make_signed_inf<float16_t>(bool sign)
{
    return u_as_f(pack_to_F16_UI(sign, 0x1F, 0));
}

template<>
inline constexpr float32_t
make_signed_inf<float32_t>(bool sign)
{
    return u_as_f(pack_to_F32_UI(sign, 0xFF, 0u));
}

template<>
inline constexpr float64_t
make_signed_inf<float64_t>(bool sign)
{
    return u_as_f(pack_to_F64_UI(sign, 0x7FF, 0u));
}

template<typename Ty>
Ty
make_signed_zero(bool sign = false);

template<>
inline constexpr uint16_t
make_signed_zero<uint16_t>(bool sign)
{
    return pack_to_F16_UI(sign, 0, 0u);
}

template<>
inline constexpr uint32_t
make_signed_zero<uint32_t>(bool sign)
{
    return pack_to_F32_UI(sign, 0, 0u);
}

template<>
inline constexpr uint64_t
make_signed_zero<uint64_t>(bool sign)
{
    return pack_to_F64_UI(sign, 0, 0u);
}

template<>
inline constexpr float16_t
make_signed_zero<float16_t>(bool sign)
{
    return u_as_f(make_signed_zero<uint16_t>(sign));
}

template<>
inline constexpr float32_t
make_signed_zero<float32_t>(bool sign)
{
    return u_as_f(make_signed_zero<uint32_t>(sign));
}

template<>
inline constexpr float64_t
make_signed_zero<float64_t>(bool sign)
{
    return u_as_f(make_signed_zero<uint64_t>(sign));
}

/**
Returns true when 32-bit unsigned integer `uiA' has the bit pattern of a
32-bit floating-point signaling NaN.
*/
inline float64_t
norm_round_pack_to_F64(bool const sign,
                       int16_t const exp,
                       uint64_t const sig)
{
    int8_t const shiftDist = count_leading_zeros(sig) - 1;
    int16_t const exp_1 = exp - shiftDist;

    if (10 <= shiftDist && static_cast<uint16_t>(exp_1) < 0x7FD) {
        return u_as_f(pack_to_F64_UI(sign, sig ? exp_1 : 0, sig << (shiftDist - 10)));
    }

    return round_pack_to_F64(sign, exp_1, sig << shiftDist);
}

inline float32_t
norm_round_pack_to_F32(bool const sign,
                       int16_t const exp,
                       uint32_t const sig)
{
    int8_t const shiftDist = count_leading_zeros(sig) - 1;
    int16_t const exp_1 = exp - shiftDist;

    if (7 <= shiftDist && static_cast<uint16_t>(exp_1) < 0xFD) {
        return u_as_f(pack_to_F32_UI(sign, sig ? exp_1 : 0, sig << (shiftDist - 7)));
    }

    return round_pack_to_F32(sign, exp_1, sig << shiftDist);
}

}  // namespace internals
}  // namespace softfloat

#endif  /* SOFTFLOAT_INTERNALS_HPP_ */
