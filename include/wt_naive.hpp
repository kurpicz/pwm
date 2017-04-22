/*******************************************************************************
 * include/wt_naive.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WT_NAIVE
#define WT_NAIVE

#include <cstring>
#include <vector>

template <typename AlphabetType, typename SizeType>
class wt_naive {

public:
  wt_naive(const std::vector<AlphabetType>& text, const SizeType size,
    const SizeType levels) : _bv(levels) {

    std::vector<AlphabetType> local_text = text;

    for (uint64_t level = 0; level < levels; ++level) {
      _bv[level] = new uint64_t[(size + 63ULL) >> 6];
      memset(_bv[level], 0, ((size + 63ULL) >> 6) * sizeof(uint64_t));

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

      std::vector<std::vector<AlphabetType>> buckets(1ULL << (level + 1));
      for (uint64_t i = 0; i < local_text.size(); ++i) {
        buckets[local_text[i] >> (levels - (level + 1))]
          .emplace_back(local_text[i]);
      }
      cur_pos = 0;
      for (const auto& bucket : buckets) {
        for (const auto character : bucket) {
          local_text[cur_pos++] = character;
        }
      }
    }
  }

  auto get_bv() const {
    return _bv;
  }

private:
  std::vector<uint64_t*> _bv;
}; // class wt_naive

#endif // WT_NAIVE

/******************************************************************************/