/*******************************************************************************
 * include/util/debug.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <iostream>
#include <string>

#include "util/wavelet_structure.hpp"

template<typename WordType, typename bv_t>
std::string bit_string(bv_t const& bv, uint64_t const size) {
  constexpr WordType BITS = (sizeof(WordType) * CHAR_BIT);

  auto s = std::string(size * BITS, '0');
  for (uint64_t bit = 0; bit < size * BITS; bit++) {
    s[bit] += bit_at<WordType>(bv, bit);
  }

  return s;
}

static std::vector<std::vector<uint64_t>> level_sizes(const base_bit_vectors& bv,
  uint64_t bit_offset, uint64_t bit_length, uint64_t level) {

  if (level == bv.levels()) {
    return {};
  }

  uint64_t zeroes = 0;

  for(uint64_t i = 0; i < bit_length; ++i) {
    if (bit_at(bv[level], bit_offset + i) == 0) {
      ++zeroes;
    }
  }

  uint64_t size_left = zeroes;
  uint64_t size_right = bit_length - zeroes;

  auto sizes_left  = level_sizes(bv, bit_offset, size_left, level + 1);
  auto sizes_right = level_sizes(
    bv, bit_offset + size_left, size_right, level + 1);

  std::vector<std::vector<uint64_t>> r;

  r.push_back({});
  r.back().push_back({bit_length});
  //r.back().push_back({size_right});

  for (uint64_t j = 0; j < sizes_left.size(); j++) {
    r.push_back({});
    auto& v = r.back();

    for(auto& e : sizes_left[j]) {
      v.push_back(e);
    }
    for(auto& e : sizes_right[j]) {
      v.push_back(e);
    }
  }

  return r;
}

/******************************************************************************/
