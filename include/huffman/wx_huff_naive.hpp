/*******************************************************************************
 * include/huffman/wx_huff_naive.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>

#include "huff_bit_vectors.hpp"
#include "huff_codes.hpp"
#include "util/wavelet_structure.hpp"

template <typename AlphabetType, bool is_tree>
class wx_huff_naive;

template <typename AlphabetType>
class wx_huff_naive<AlphabetType, true> {

public:
  static constexpr bool    is_parallel = false;
  static constexpr bool    is_tree     = true;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);
  static constexpr bool    is_huffman_shaped = true;

  static wavelet_structure compute(AlphabetType const* const text,
    const uint64_t size, const uint64_t /*levels*/) {

    canonical_huff_codes<AlphabetType, is_tree> codes(text, size);

    if(size == 0) {
      return wavelet_structure_tree_huffman<AlphabetType>(std::move(codes));
    }

    const uint64_t levels = codes.levels();

    auto _bv = huff_bit_vectors(levels, codes.level_sizes());
    auto& bv = _bv.raw_data();

    std::vector<AlphabetType> local_text(size);
    for(size_t i = 0; i < size; i++) {
      local_text[i] = text[i];
    }

    for (uint64_t level = 0; level < levels; ++level) {
      uint32_t cur_pos = 0;
      for (; cur_pos + 64 <= local_text.size(); cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < 64; ++i) {
          word <<= 1;
          word |= codes.encode_symbol(local_text[cur_pos + i])[level];
        }
        bv[level][cur_pos >> 6] = word;
      }
      if (local_text.size() & 63ULL) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < local_text.size() - cur_pos; ++i) {
          word <<= 1;
          word |= codes.encode_symbol(local_text[cur_pos + i])[level];
        }
        word <<= (64 - (local_text.size() & 63ULL));
        bv[level][local_text.size() >> 6] = word;
      }
      std::vector<std::vector<AlphabetType>> buckets(1ULL << (level + 1));
      for (uint64_t i = 0; i < local_text.size(); ++i) {
        const AlphabetType cur_symbol = local_text[i];
        if (codes.code_length(cur_symbol) > level + 1) {
          buckets[codes.encode_symbol(cur_symbol).prefix(level + 1)]
            .emplace_back(cur_symbol);
        }
      }
      for (uint64_t i = 1; i < buckets.size(); ++i) {
        std::move(buckets[i].begin(), buckets[i].end(),
          std::back_inserter(buckets[0]));
      }
      local_text.swap(buckets[0]);
    }
    return wavelet_structure_tree_huffman<AlphabetType>(std::move(_bv), std::move(codes));
  }
}; // class wt_huff_naive

template <typename AlphabetType>
class wx_huff_naive<AlphabetType, false> {

public:
  static constexpr bool    is_parallel = false;
  static constexpr bool    is_tree     = false;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);
  static constexpr bool    is_huffman_shaped = true;

  static wavelet_structure compute(AlphabetType const* const text,
    const uint64_t size, const uint64_t /*levels*/) {

    canonical_huff_codes<AlphabetType, is_tree> codes(text, size);

    if(size == 0) {
      return wavelet_structure_matrix_huffman<AlphabetType>(std::move(codes));
    }

    const uint64_t levels = codes.levels();

    auto _bv = huff_bit_vectors(levels, codes.level_sizes());
    auto& bv = _bv.raw_data();
    auto _zeros = std::vector<size_t>(levels, 0);

    std::vector<AlphabetType> local_text(size);
    for(size_t i = 0; i < size; i++) {
      local_text[i] = text[i];
    }

    // Construct each level top-down
    for (uint64_t level = 0; level < levels; ++level) {
      // Insert the level-th MSB in the bit vector of the level (in text order)
      uint32_t cur_pos = 0;
      for (; cur_pos + 64 <= local_text.size(); cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < 64; ++i) {
          word <<= 1;
          word |= codes.encode_symbol(local_text[cur_pos + i])[level];
        }
        bv[level][cur_pos >> 6] = word;
      }
      if (local_text.size() & 63ULL) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < local_text.size() - cur_pos; ++i) {
          word <<= 1;
          word |= codes.encode_symbol(local_text[cur_pos + i])[level];
        }
        word <<= (64 - (local_text.size() & 63ULL));
        bv[level][local_text.size() >> 6] = word;
      }

      std::vector<AlphabetType> text0;
      std::vector<AlphabetType> text1;
      // Scan the text and separate characters that inserted 0s and 1s
      for (const auto symbol : local_text) {
        const code_pair cp = codes.encode_symbol(symbol);
        if (cp.code_length > level + 1) {
          if (cp[level]) {
            text1.emplace_back(symbol);
          } else {
            text0.emplace_back(symbol);
          }
        }
      }
      _zeros[level] = text0.size();
      std::move(text1.begin(), text1.end(), std::back_inserter(text0));
      local_text.swap(text0);
    }
    return wavelet_structure_matrix_huffman<AlphabetType>(std::move(_bv), std::move(_zeros), std::move(codes));
  }
}; // class wx_huff_naive<MATRIX>

/******************************************************************************/
