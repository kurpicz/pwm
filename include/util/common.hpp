/*******************************************************************************
 * include/common.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef COMMON
#define COMMON

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

inline auto rho_identity(uint64_t levels) {
    return [](auto level, auto i) -> uint64_t {
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
    static auto create(size_t levels) {
        return rho_bit_reverse(levels);
    }
};

template<>
struct rho_dispatch<false> {
    static auto create(size_t levels) {
        return rho_identity(levels);
    }
};

constexpr size_t word_size(uint64_t size) {
    return (size + 63ULL) >> 6;
}

class Bvs {
    std::vector<uint64_t*> m_data;
    uint64_t m_size;
public:
    inline Bvs(): m_size(0) {}
    inline Bvs(uint64_t size, uint64_t levels):
        m_data(levels), m_size(size)
    {
        assert(levels != 0);

        m_data[0] = new uint64_t[word_size(size) * levels];
        memset(m_data[0], 0, (word_size(size) * sizeof(uint64_t)) * levels);

        for (uint64_t level = 1; level < levels; ++level) {
            m_data[level] = m_data[level - 1] + word_size(size);
        }
    }

    inline uint64_t levels() const {
        return m_data.size();
    }

    inline uint64_t size() const {
        return m_size;
    }

    inline const std::vector<uint64_t*>& vec() const {
        return m_data;
    }

    inline std::vector<uint64_t*>& vec() {
        return m_data;
    }

    inline ~Bvs() {
        if (m_data.size() > 0) {
            delete[] m_data[0];
        }
    }

    inline Bvs(Bvs&& other):
        m_data(std::move(other.m_data)),
        m_size(other.m_size) {}

    inline Bvs& operator=(Bvs&& other) {
        if (m_data.size() > 0) {
            delete[] m_data[0];
        }
        m_data = std::move(other.m_data);
        m_size = other.m_size;

        return *this;
    }

};

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

#endif // COMMON

/******************************************************************************/
