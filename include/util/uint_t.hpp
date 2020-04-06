/*******************************************************************************
 * include/uint_t.hpp
 *
 * Copyright (C) 2016 Marvin Löbel
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/uint_t/IntegerBase.hpp"

#include <cstdint>
#include <cmath>

namespace pwm {

template<size_t bits>
class uint_impl_t;

template<class MB>
struct UinttDispatch {
    typedef MB SelfMaxBit;

    template<class Ref, class V>
    inline static void assign(Ref& self, V v) {
        self.m_data = v;
    }

    template<class Ref, class R>
    inline static R cast_for_op(const Ref& self) {
        return self.m_data;
    }
};

template<size_t N>
struct ConstIntegerBaseTrait<uint_impl_t<N>, typename std::enable_if<(N <= 32)>::type> {
    typedef UinttDispatch<uint32_t> Dispatch;
};

template<size_t N>
struct IntegerBaseTrait<uint_impl_t<N>, typename std::enable_if<(N <= 32)>::type>
: ConstIntegerBaseTrait<uint_impl_t<N>> {
    typedef UinttDispatch<uint32_t> Dispatch;
};

template<size_t N>
struct ConstIntegerBaseTrait<uint_impl_t<N>, typename std::enable_if<(N > 32)>::type> {
    typedef UinttDispatch<uint64_t> Dispatch;
};

template<size_t N>
struct IntegerBaseTrait<uint_impl_t<N>, typename std::enable_if<(N > 32)>::type>
: ConstIntegerBaseTrait<uint_impl_t<N>> {
    typedef UinttDispatch<uint64_t> Dispatch;
};

template<>
struct ConstIntegerBaseTrait<bool> {
    typedef UinttDispatch<uint32_t> Dispatch;
};

template<>
struct IntegerBaseTrait<bool>
: ConstIntegerBaseTrait<bool> {
    typedef UinttDispatch<uint32_t> Dispatch;
};

/// Custom integer type for storing values of arbitrary bit size `bits`.
///
/// It is guaranteed that this type will only have a byte size as large
/// as needed to store all bits, but there will be padding bits for values
/// not multiples of 8.
///
/// In practice useful sizes are 40, 48 and 56 bits.
/// 40 bit indices correspond to 1TiB of addressable memory,
/// and 48 bit correspond to today's hardware limits of the x86_64 architecture.
template<size_t bits>
class uint_impl_t: public IntegerBase<uint_impl_t<bits>> {
    static_assert(bits > 0, "bits must be non-negative");
    static_assert(bits < 65, "bits must be at most 64");
    uint64_t m_data: bits;

    friend struct UinttDispatch<uint32_t>;
    friend struct UinttDispatch<uint64_t>;

public:
    constexpr uint_impl_t(): m_data(0) {}
    constexpr uint_impl_t(uint_impl_t&& i): m_data(i.m_data) {}

    // copying
    constexpr uint_impl_t(const uint_impl_t& i): m_data(i.m_data) {}
    inline uint_impl_t& operator=(const uint_impl_t& b) { m_data = b.m_data; return *this; }

    // conversions for all fundamental char types
    constexpr uint_impl_t(unsigned char i): m_data(i) {}
    inline uint_impl_t& operator=(unsigned char data) { m_data = data; return *this; }
    constexpr operator unsigned char() const { return m_data; }

    constexpr uint_impl_t(signed char i): m_data(i) {}
    inline uint_impl_t& operator=(signed char data) { m_data = data; return *this; }
    constexpr operator signed char() const { return m_data; }

    constexpr uint_impl_t(char i): m_data(i) {}
    constexpr uint_impl_t& operator=(char data) { m_data = data; return *this; }
    constexpr operator char() const { return m_data; }

    // conversions for all fundamental integer types
    constexpr uint_impl_t(unsigned int i): m_data(i) {}
    inline uint_impl_t& operator=(unsigned int data) { m_data = data; return *this; }
    constexpr operator unsigned int() const { return m_data; }

    constexpr uint_impl_t(unsigned short int i): m_data(i) {}
    inline uint_impl_t& operator=(unsigned short int data) { m_data = data; return *this; }
    constexpr operator unsigned short int() const { return m_data; }

    constexpr uint_impl_t(unsigned long int i): m_data(i) {}
    inline uint_impl_t& operator=(unsigned long int data) { m_data = data; return *this; }
    constexpr operator unsigned long int() const { return m_data; }

    constexpr uint_impl_t(unsigned long long int i): m_data(i) {}
    inline uint_impl_t& operator=(unsigned long long int data) { m_data = data; return *this; }
    constexpr operator unsigned long long int() const { return m_data; }

    constexpr uint_impl_t(int i): m_data(i) {}
    inline uint_impl_t& operator=(int data) { m_data = data; return *this; }
    constexpr operator int() const { return m_data; }

    constexpr uint_impl_t(short int i): m_data(i) {}
    inline uint_impl_t& operator=(short int data) { m_data = data; return *this; }
    constexpr operator short int() const { return m_data; }

    constexpr uint_impl_t(long int i): m_data(i) {}
    inline uint_impl_t& operator=(long int data) { m_data = data; return *this; }
    constexpr operator long int() const { return m_data; }

    constexpr uint_impl_t(long long int i): m_data(i) {}
    inline uint_impl_t& operator=(long long int data) { m_data = data; return *this; }
    constexpr operator long long int() const { return m_data; }
} __attribute__((packed));

template<size_t N>
inline std::ostream& operator<<(std::ostream& os, const uint_impl_t<N>& v) {
    return os << uint64_t(v);
}

template<size_t N>
inline std::istream& operator>>(std::istream& is, uint_impl_t<N>& v) {
    uint64_t v2;
    auto& x = is >> v2;
    v = v2;
    return x;
}

// Specialize uint_t<1> to be identical to bool

// We disable 1-Bit instantiation of the actual type because the typedef
// below will redirect uint_t<1> to bool.
template<>
class uint_impl_t<1> {};

template<size_t N>
struct uint_dispatch_t {
    using type = uint_impl_t<N>;
};

template<>
struct uint_dispatch_t<1> {
    using type = bool;
};

template<size_t N>
using uint_t = typename uint_dispatch_t<N>::type;

static_assert(sizeof(uint_t<8>)  == 1, "sanity check");
static_assert(sizeof(uint_t<16>) == 2, "sanity check");
static_assert(sizeof(uint_t<24>) == 3, "sanity check");
static_assert(sizeof(uint_t<32>) == 4, "sanity check");
static_assert(sizeof(uint_t<40>) == 5, "sanity check");
static_assert(sizeof(uint_t<48>) == 6, "sanity check");
static_assert(sizeof(uint_t<56>) == 7, "sanity check");
static_assert(sizeof(uint_t<64>) == 8, "sanity check");

static_assert(sizeof(uint_t<7>)  == 1, "sanity check");
static_assert(sizeof(uint_t<15>) == 2, "sanity check");
static_assert(sizeof(uint_t<23>) == 3, "sanity check");
static_assert(sizeof(uint_t<31>) == 4, "sanity check");
static_assert(sizeof(uint_t<39>) == 5, "sanity check");
static_assert(sizeof(uint_t<47>) == 6, "sanity check");
static_assert(sizeof(uint_t<55>) == 7, "sanity check");
static_assert(sizeof(uint_t<63>) == 8, "sanity check");

static_assert(sizeof(uint_t<9>)  == 2, "sanity check");
static_assert(sizeof(uint_t<17>) == 3, "sanity check");
static_assert(sizeof(uint_t<25>) == 4, "sanity check");
static_assert(sizeof(uint_t<33>) == 5, "sanity check");
static_assert(sizeof(uint_t<41>) == 6, "sanity check");
static_assert(sizeof(uint_t<49>) == 7, "sanity check");
static_assert(sizeof(uint_t<57>) == 8, "sanity check");

}

using uint24_t = pwm::uint_t<24>;
using uint40_t = pwm::uint_t<40>;
using uint48_t = pwm::uint_t<48>;
using uint56_t = pwm::uint_t<56>;

namespace std {
    template<size_t N>
    class numeric_limits<pwm::uint_impl_t<N>> {
        using T = pwm::uint_impl_t<N>;
    public:
        static constexpr bool is_specialized = true;
        static constexpr bool is_signed = false;
        static constexpr bool is_integer = true;
        static constexpr bool is_exact = true;
        static constexpr bool has_infinity = false;
        static constexpr bool has_quiet_NaN = false;
        static constexpr bool has_signaling_NaN = false;
        static constexpr std::float_denorm_style has_denorm = std::denorm_absent;
        static constexpr bool has_denorm_loss = false;
        static constexpr std::float_round_style round_style = std::round_toward_zero;
        static constexpr bool is_iec559 = false;
        static constexpr bool is_bounded = true;

        // TODO: Not true currently, would need code to support it
        static constexpr bool is_modulo = false;

        static constexpr int digits = N;
        static constexpr int digits10 = std::numeric_limits<T>::digits * std::log10(2);
        static constexpr int max_digits10 = 0;
        static constexpr int radix = 2;
        static constexpr int min_exponent = 0;
        static constexpr int min_exponent10 = 0;
        static constexpr int max_exponent = 0;
        static constexpr int max_exponent10 = 0;

        // Maybe need to derive from base uint64_t type
        static constexpr bool traps = true;
        static constexpr bool tinyness_before = false;

        static constexpr T min() { return 0; }
        static constexpr T lowest() { return 0; }
        static constexpr T max() { return uint64_t(std::pow(2, N) - 1); }
        static constexpr T epsilon() { return 0; }
        static constexpr T round_error() { return 0; }
        static constexpr T infinity() { return 0; }
        static constexpr T quiet_NaN() { return 0; }
        static constexpr T signaling_NaN() { return 0; }
        static constexpr T denorm_min() { return 0; }
    };

    template<size_t N>
    std::string to_string(pwm::uint_impl_t<N> value) {
        return std::to_string(uint64_t(value));
    }

    template<size_t N>
    struct hash<pwm::uint_impl_t<N>> {
        size_t operator()(const pwm::uint_impl_t<N>& x) const {
            return hash<uint64_t>()(x);
        }
    };
}
