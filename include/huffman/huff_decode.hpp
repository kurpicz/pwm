/*******************************************************************************
 * include/huffman/huff_decode.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "huff_bit_vectors.hpp"
#include "util/debug.hpp"

[[gnu::unused]] // TODO: C++17 [[maybe_unused]] 
static std::string decode_wt_huff(const huff_bit_vectors& bv) {
  std::vector<uint64_t> code_words;
  std::vector<uint64_t> code_lengths;

  std::vector<std::vector<uint64_t>> borders(bv.levels());

  //  Count initial 0s and 1s.
  for (uint64_t i = 0; i < bv.level_size(0); ++i) {

  }

  for (uint64_t i = 0; i + 1 < bv.levels(); ++i) {
    // for (uint64_t j = 0; j <)
  } 

  std::vector<uint8_t> r;
  return std::string(r.begin(), r.end());
}

[[gnu::unused]] // TODO: C++17 [[maybe_unused]] 
static std::string decode_wm_huff(const huff_bit_vectors& bv,
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
