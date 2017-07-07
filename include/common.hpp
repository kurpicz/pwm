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

#endif // COMMON

/******************************************************************************/
