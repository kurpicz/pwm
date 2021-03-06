/*******************************************************************************
 * include/wx_naive.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>

#include "arrays/bit_vectors.hpp"
#include "construction/wavelet_structure.hpp"

#include "wx_base.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_naive;

template <typename AlphabetType>
class wx_naive<AlphabetType, true> : public wx_in_out_external<false, false> {

public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree = true;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  static wavelet_structure compute(AlphabetType const* const text,
                                   const uint64_t size,
                                   const uint64_t levels) {

    if (size == 0) {
      return wavelet_structure_tree();
    }

    auto bv = bit_vectors<>(levels, size);

    // TODO: When not semi_external, this can/should be done way nicer.
    std::vector<AlphabetType> local_text(size);
    memcpy(local_text.data(), text, size * sizeof(AlphabetType));

    for (uint64_t level = 0; level < levels; ++level) {
      uint32_t cur_pos = 0;
      for (; cur_pos + 64 <= size; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < 64; ++i) {
          word <<= 1;
          word |= ((local_text[cur_pos + i] >> (levels - (level + 1))) & 1ULL);
        }
        bv[level][cur_pos >> 6] = word;
      }
      if (size & 63ULL) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < size - cur_pos; ++i) {
          word <<= 1;
          word |= ((local_text[cur_pos + i] >> (levels - (level + 1))) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        bv[level][size >> 6] = word;
      }
      if (level + 1 < levels) {
        std::vector<std::vector<AlphabetType>> buckets(1ULL << (level + 1));
        for (uint64_t i = 0; i < local_text.size(); ++i) {
          buckets[local_text[i] >> (levels - (level + 1))].emplace_back(
              local_text[i]);
        }
        cur_pos = 0;
        for (const auto& bucket : buckets) {
          for (const auto character : bucket) {
            local_text[cur_pos++] = character;
          }
        }
      }
    }
    return wavelet_structure_tree(std::move(bv));
  }
}; // class wt_naive

template <typename AlphabetType>
class wx_naive<AlphabetType, false> : public wx_in_out_external<false, false> {

public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree = false;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  static wavelet_structure compute(AlphabetType const* const text,
                                   const uint64_t size,
                                   const uint64_t levels) {

    if (size == 0) {
      return wavelet_structure_matrix();
    }

    auto bv = bit_vectors<>(levels, size);
    auto _zeros = std::vector<size_t>(levels, 0);

    // TODO: When not semi_xternal, this can/should be done way nicer.
    std::vector<AlphabetType> local_text(size);
    for (size_t i = 0; i < size; i++) {
      local_text[i] = text[i];
    }

    // Construct each level top-down
    for (uint64_t level = 0; level < levels; ++level) {
      // Insert the level-th MSB in the bit vector of the level (in text order)
      uint32_t cur_pos = 0;
      for (; cur_pos + 64 <= size; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < 64; ++i) {
          word <<= 1;
          word |= ((local_text[cur_pos + i] >> (levels - (level + 1))) & 1ULL);
        }
        bv[level][cur_pos >> 6] = word;
      }
      if (size & 63ULL) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < size - cur_pos; ++i) {
          word <<= 1;
          word |= ((local_text[cur_pos + i] >> (levels - (level + 1))) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        bv[level][size >> 6] = word;
      }

      std::vector<AlphabetType> text0;
      text0.reserve(size);
      std::vector<AlphabetType> text1;
      text1.reserve(size);
      // Scan the text and separate characters that inserted 0s and 1s
      for (uint64_t i = 0; i < size; ++i) {
        if ((local_text[i] >> (levels - (level + 1))) & 1ULL) {
          text1.push_back(local_text[i]);
        } else {
          text0.push_back(local_text[i]);
        }
      }
      _zeros[level] = text0.size();
      if (level + 1 < levels) {
        // "Sort" the text stably based on the bit inserted in the bit vector
        for (uint64_t i = 0; i < text0.size(); ++i) {
          local_text[i] = text0[i];
        }
        for (uint64_t i = 0; i < text1.size(); ++i) {
          local_text[i + text0.size()] = text1[i];
        }
      }
    }
    return wavelet_structure_matrix(std::move(bv), std::move(_zeros));
  }
}; // class wx_naive<MATRIX>

/******************************************************************************/
