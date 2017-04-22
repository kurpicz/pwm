/*******************************************************************************
 * include/wm_naive.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_NAIVE
#define WM_NAIVE

#include <cstring>
#include <vector>

#include "common.hpp"

template <typename AlphabetType, typename SizeType>
class wm_naive {

public:
  wm_naive(const std::vector<AlphabetType>& text, const SizeType size,
    const SizeType levels) : _bv(levels), _zeros(levels, 0) {

    std::vector<AlphabetType> local_text = text;

    // Construct each level top-down
    for (uint64_t level = 0; level < levels; ++level) {
      _bv[level] = new uint64_t[(size + 63ULL) >> 6];
      memset(_bv[level], 0, ((size + 63ULL) >> 6) * sizeof(uint64_t));

      // Insert the level-th MSB in the bit vector of the level (in text order)
      uint32_t cur_pos = 0;
      for (; cur_pos + 64 <= size; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < 64; ++i) {
          word <<= 1;
          word |= ((local_text[cur_pos + i] >> (levels - (level + 1))) & 1ULL);
        }
        _bv[level][cur_pos >> 6] = word;
      }
      if (size & 63ULL) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < size - cur_pos; ++i) {
          word <<= 1;
          word |= ((local_text[cur_pos + i] >> (levels - (level + 1))) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        _bv[level][size >> 6] = word;
      }

      std::vector<AlphabetType> text0;
      text0.reserve(size);
      std::vector<AlphabetType> text1;
      text1.reserve(size);
      // Scan the text and separate characters that inserted 0s and 1s
      for (uint64_t i = 0; i < local_text.size(); ++i) {
        if ((local_text[i] >> (levels - (level + 1))) & 1ULL) {
          text1.push_back(local_text[i]);
        } else {
          text0.push_back(local_text[i]);
        }
      }
      // "Sort" the text stably based on the bit inserted in the bit vector
      for (uint64_t i = 0; i < text0.size(); ++i) {
        local_text[i] = text0[i];
      }
      for (uint64_t i = 0; i < text1.size(); ++i) {
        local_text[i + text0.size()] = text1[i];
      }
      _zeros[level] = text0.size();
    }
  }

  auto get_bv_and_zeros() const {
    return std::make_pair(_bv, _zeros);
  }

private:
  std::vector<uint64_t*> _bv;
  std::vector<SizeType> _zeros;
}; // class wm_naive

#endif // WM_NAIVE

/******************************************************************************/