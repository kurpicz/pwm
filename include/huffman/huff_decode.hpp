/*******************************************************************************
 * include/huffman/huff_decode.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>

#include "bit_vector/bit_vectors.hpp"
#include "huffman/huff_codes.hpp"
#include "util/border_hist_array.hpp"
#include "util/debug.hpp"
#include "util/wavelet_structure.hpp"

template<typename AlphabetType>
std::string decode_wt_huff(const base_bit_vectors& bv,
  const canonical_huff_codes<AlphabetType, true>& codes) {

  if (bv.levels() == 0) {
    return { };
  }

  struct icp {
    uint64_t index;
    code_pair cp;
  }; // struct icp

  // All codes have length greater or equal to 1. By choosing the initial length
  // to be 0, we can use the length to mark finished codes. To this end, we
  // interpret every code_length > 0 as the code_length of a finished code.
  std::vector<icp> icps(bv.level_bit_size(0), { 0, { 0ULL, 0ULL } });
  for (uint64_t i = 0; i < icps.size(); ++i) {
    icps[i].index = i;
  }

  for (uint64_t level = 0; level < bv.levels(); ++level) {
    uint64_t bit_pos = 0;
    for (; bit_pos < bv.level_bit_size(level); ++bit_pos) {
      icps[bit_pos].cp.code_word <<= 1;
      icps[bit_pos].cp.code_word |= bit_at(bv[level], bit_pos);
    }
    while (bit_pos < bv.level_bit_size(((level == 0) ? 0 : level - 1))) {
      icps[bit_pos++].cp.code_length = level;
    }
    std::stable_sort(icps.begin(), icps.begin() + bv.level_bit_size(level),
      [](const icp& a, const icp& b) {
        return a.cp.code_word < b.cp.code_word;
      });
  }
  for (uint64_t i = 0; i < bv.level_bit_size(bv.levels() - 1); ++i) {
    icps[i].cp.code_length = bv.levels();
  }

  std::sort(icps.begin(), icps.end(), [](const icp& a, const icp& b) {
    return a.index < b.index;
  });

  std::string result;
  for (const auto& icp : icps) {
    result.push_back(codes.decode_symbol(icp.cp.code_length, icp.cp.code_word));
  }
  return result;
}

template<typename AlphabetType>
std::string decode_wm_huff(const base_bit_vectors& bv,
  const canonical_huff_codes<AlphabetType, false>& codes) {

  if (bv.levels() == 0) {
    return { };
  }

  struct icp {
    uint64_t index;
    code_pair cp;
  }; // struct icp

  // All codes have length greater or equal to 1. By choosing the initial length
  // to be 0, we can use the length to mark finished codes. To this end, we
  // interpret every code_length > 0 as the code_length of a finished code.
  std::vector<icp> icps(bv.level_bit_size(0), { 0, { 0ULL, 0ULL } });
  for (uint64_t i = 0; i < icps.size(); ++i) {
    icps[i].index = i;
  }

  for (uint64_t level = 0; level < bv.levels(); ++level) {
    uint64_t bit_pos = 0;
    for (; bit_pos < bv.level_bit_size(level); ++bit_pos) {
      icps[bit_pos].cp.code_word <<= 1;
      icps[bit_pos].cp.code_word |= bit_at(bv[level], bit_pos);
    }
    while (bit_pos < bv.level_bit_size(((level == 0) ? 0 : level - 1))) {
      icps[bit_pos++].cp.code_length = level;
    }
    std::stable_sort(icps.begin(), icps.begin() + bv.level_bit_size(level),
      [](const icp& a, const icp& b) {
        return (a.cp.code_word & 1ULL) < (b.cp.code_word & 1ULL);
      });
  }
  for (uint64_t i = 0; i < bv.level_bit_size(bv.levels() - 1); ++i) {
    icps[i].cp.code_length = bv.levels();
  }

  std::sort(icps.begin(), icps.end(), [](const icp& a, const icp& b) {
    return a.index < b.index;
  });

  std::string result;
  for (const auto& icp : icps) {
    result.push_back(codes.decode_symbol(icp.cp.code_length, icp.cp.code_word));
  }
  return result;
}

/******************************************************************************/
