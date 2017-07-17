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

#include <vector>
#include <cassert>
#include <type_traits>
#include <climits>

template <typename SizeType>
static inline std::vector<SizeType> BitReverse(const SizeType levels) {
  std::vector<SizeType> result(1 << levels);
  result[0] = 0;
  result[1] = 1;
  for (SizeType i = 1; i < levels; ++i) {
    for (SizeType j = 0; j < (1u << i); ++j) {
      result[j] <<= 1;
    }
    for (SizeType j = 0; j < (1u << i); ++j) {
      result[j + (1 << i)] = result[j] + 1;
    }
  }
  return result;
}

template <typename SizeType>
static inline void BitReverse(const SizeType levels, SizeType* out) {
  out[0] = 0;
  out[1] = 1;
  for (SizeType i = 1; i < levels; ++i) {
    for (SizeType j = 0; j < (1 << i); ++j) {
      out[j] <<= 1;
    }
    for (SizeType j = 0; j < (1 << i); ++j) {
      out[j + (1 << i)] = out[j] + 1;
    }
  }
}

template <typename SizeType>
inline auto rho_identity(SizeType levels) {
    return [](auto level, auto i) -> SizeType {
        return i;
    };
}

template <typename SizeType>
inline auto rho_bit_reverse(SizeType levels) {
    auto bit_reverse = std::vector<std::vector<SizeType>>(levels);
    bit_reverse[levels - 1] = BitReverse<SizeType>(levels - 1);
    for(size_t level = levels - 1; level > 0; level--) {
        bit_reverse[level - 1] = std::vector<SizeType>(bit_reverse[level].size() / 2);
        for(size_t i = 0; i < bit_reverse[level - 1].size(); i++) {
            bit_reverse[level - 1][i] = bit_reverse[level][i] >> 1;
        }
    }

    return [bit_reverse = std::move(bit_reverse)](auto level, auto i) -> SizeType {
        return bit_reverse[level][i];
    };
}

constexpr size_t word_size(uint64_t size) {
    return (size + 63ULL) >> 6;
}

template<typename SizeType>
class Bvs {
    std::vector<uint64_t*> m_data;
    SizeType m_size;
public:
    inline Bvs(): m_size(0) {}
    inline Bvs(SizeType size, SizeType levels):
        m_data(levels)
    {
        assert(levels != 0);

        m_data[0] = new uint64_t[word_size(size) * levels];
        memset(m_data[0], 0, (word_size(size) * sizeof(uint64_t)) * levels);

        for (SizeType level = 1; level < levels; ++level) {
            m_data[level] = m_data[level - 1] + word_size(size);
        }
    }

    inline SizeType levels() const {
        return m_data.size();
    }

    inline SizeType size() const {
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
        m_size(std::move(other.m_size)) {}

    inline Bvs& operator=(Bvs&& other) {
        if (m_data.size() > 0) {
            delete[] m_data[0];
        }
        m_data = std::move(other.m_data);
        m_size = std::move(other.m_size);

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
