/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <stdint.h>
#include <vector>
#include <cassert>
#include <type_traits>
#include <climits>
#include <cstring>

using permutation_type = std::vector<uint64_t> (*)(const uint64_t levels);

static inline std::vector<uint64_t> BitReverse(const uint64_t levels) {
  std::vector<uint64_t> result(1 << levels);
  result[0] = 0;
  result[1] = 1;
  for (uint64_t i = 1; i < levels; ++i) {
    for (uint64_t j = 0; j < (1u << i); ++j) {
      result[j] <<= 1;
    }
    for (uint64_t j = 0; j < (1u << i); ++j) {
      result[j + (1 << i)] = result[j] + 1;
    }
  }
  return result;
}

static inline void BitReverse(const uint64_t levels, uint64_t* out) {
  out[0] = 0;
  out[1] = 1;
  for (uint64_t i = 1; i < levels; ++i) {
    for (uint64_t j = 0; j < (1ULL << i); ++j) {
      out[j] <<= 1;
    }
    for (uint64_t j = 0; j < (1ULL << i); ++j) {
      out[j + (1 << i)] = out[j] + 1;
    }
  }
}

inline auto rho_identity(uint64_t /*levels*/) {
    return [](auto /*level*/, auto i) -> uint64_t {
        return i;
    };
}

inline auto rho_bit_reverse(uint64_t levels) {
    auto bit_reverse = std::vector<std::vector<uint64_t>>(levels);
    bit_reverse[levels - 1] = BitReverse(levels - 1);
    for(size_t level = levels - 1; level > 0; level--) {
        bit_reverse[level - 1] = std::vector<uint64_t>(bit_reverse[level].size() / 2);
        for(size_t i = 0; i < bit_reverse[level - 1].size(); i++) {
            bit_reverse[level - 1][i] = bit_reverse[level][i] >> 1;
        }
    }

    return [bit_reverse = std::move(bit_reverse)](auto level, auto i) -> uint64_t {
        return bit_reverse[level][i];
    };
}

template<bool is_matrix>
struct rho_dispatch {};

template<>
struct rho_dispatch<true> {
    using type = decltype(rho_bit_reverse(0));

    static auto create(size_t levels) {
        return rho_bit_reverse(levels);
    }
};

template<>
struct rho_dispatch<false> {
    using type = decltype(rho_identity(0));

    static auto create(size_t levels) {
        return rho_identity(levels);
    }
};

constexpr uint64_t word_size(uint64_t size) {
  return (size + 63ULL) >> 6;
}

// A little template helper for dropping a type early
template <typename T>
void drop_me(T const&) = delete;
template <typename T>
void drop_me(T&) = delete;
template <typename T>
void drop_me(T const&&) = delete;
template <typename T>
void drop_me(T&& t) {
    std::remove_reference_t<T>(std::move(t));
}

constexpr size_t log2(size_t n) {
    return (n < 2) ? 1 : 1 + log2(n / 2);
}

template<typename WordType = uint64_t, typename bv_t>
inline auto bit_at(const bv_t& bv, size_t i) -> bool {
    constexpr WordType BITS = (sizeof(WordType) * CHAR_BIT);
    constexpr WordType MOD_MASK = BITS - 1;

    size_t offset = i / BITS;
    size_t word_offset = i & MOD_MASK;
    return (bv[offset] >> (MOD_MASK - word_offset)) & 1ull;
}

/******************************************************************************/
