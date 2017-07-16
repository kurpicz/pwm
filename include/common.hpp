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
        DCHECK_GT(levels, 0);

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

#endif // COMMON

/******************************************************************************/
