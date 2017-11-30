/*******************************************************************************
 * include/huffman/huff_decode.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/debug.hpp"

[[gnu::unused]] // TODO: C++17 [[maybe_unused]] 
static std::string decode_wt_huff(const bit_vectors& bv, uint64_t length) {
  auto ls = level_sizes(bv, 0, length, 0);

  for (auto& v : ls) {
    for (uint64_t i = 1; i < v.size(); i++) {
      v[i] = v[i - 1] + v[i];
    }
    for (uint64_t i = 1; i < v.size(); i++) {
      uint64_t j = v.size() - i;
      v[j] = v[j - 1];
    }
    if (v.size() > 0) {
      v[0] = 0;
    }
  }

  auto r = std::vector<uint8_t>(length);

  for (uint64_t i = 0; i < length; i++) {
    uint8_t value = 0;
    uint64_t j = 0;
    for (uint64_t level = 0; level < bv.levels(); level++) {
      auto& offset = ls[level][j];
      uint8_t bit = bit_at(bv[level], offset);

      value <<= 1;
      value |= bit;

      offset++;
      j = 2 * j + bit;
    }
    r[i] = value;
  }

  return std::string(r.begin(), r.end());
}

[[gnu::unused]] // TODO: C++17 [[maybe_unused]] 
static std::string decode_wm_huff(const bit_vectors& bv,
  const std::vector<uint64_t>& zeros, const uint64_t length) {

  if (bv.levels() == 0) {
    return {};
  }

  auto r = std::vector<uint8_t>(length, uint8_t(0));
  auto rtmp = std::vector<uint8_t>(length, uint8_t(0));
  // print_bv_zeros(bv, zeros, length);

  for(uint64_t level = bv.levels() - 1; level > 0; level--) {
    uint64_t offset0 = 0;
    uint64_t offset1 = zeros[level - 1];

    for(uint64_t i = 0; i < length; i++) {
      r[i] |= (bit_at(bv[level], i) << (bv.levels() - level - 1));
    }

    for(uint64_t i = 0; i < length; i++) {
      if(bit_at(bv[level - 1], i) == 0) {
        rtmp[i] = r[offset0];
        offset0++;
      } else {
        rtmp[i] = r[offset1];
        offset1++;
      }
    }

    r.swap(rtmp);
  }

  for(uint64_t i = 0; i < length; i++) {
    r[i] |= bit_at(bv[0], i) << (bv.levels() - 1);
  }

  return std::string(r.begin(), r.end());
}

/******************************************************************************/
