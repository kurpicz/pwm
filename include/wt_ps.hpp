/*******************************************************************************
 * include/wt_ps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WT_PREFIX_SORTING
#define WT_PREFIX_SORTING

#include <vector>

#include "common.hpp"

template <typename AlphabetType, typename SizeType>
class wt_ps {

public:
  wt_ps(const std::vector<AlphabetType>& text, const SizeType size,
    const SizeType levels) : _bv(levels) {

    SizeType cur_max_char = (1 << levels);
    std::vector<SizeType> s_pos(cur_max_char, 0);
    std::vector<SizeType> hist(cur_max_char, 0);
    std::vector<SizeType> borders(cur_max_char, 0);
    std::vector<AlphabetType> sorted_text(size);

    _bv[0] = new uint64_t[(size + 63ULL) >> 6];
    memset(_bv[0], 0, ((size + 63ULL) >> 3)); // memset is ok (all to 0)

    SizeType cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (SizeType i = 0; i < 64; ++i) {
        ++hist[text[cur_pos + i]];
        word <<= 1;
        word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      }
      (_bv[0])[cur_pos >> 6] = word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (SizeType i = 0; i < size - cur_pos; ++i) {
        ++hist[text[cur_pos + i]];
        word <<= 1;
        word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      }
      word <<= (64 - (size & 63ULL));
      (_bv[0])[size >> 6] = word;
    }

    for (SizeType level = levels - 1; level > 0; --level) {
      _bv[level] = new uint64_t[(size + 63ULL) >> 6];
      memset(_bv[level], 0, ((size + 63ULL) >> 3)); // memset is ok (all to 0)

      cur_max_char >>= 1;
      for (SizeType i = 0; i < cur_max_char; ++i) {
        hist[i] = hist[i << 1] + hist[(i << 1) + 1];
      }
      borders[0] = 0;
      for (SizeType i = 1; i < cur_max_char; ++i) {
        borders[i] = borders[i - 1] + hist[i - 1];
      }

      for (SizeType i = 0; i < size; ++i) {
        const AlphabetType cur_char = text[i];
        sorted_text[borders[cur_char >> (levels - level)]++] = cur_char;
      }

      cur_pos = 0;
      for (; cur_pos + 63 < size; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (SizeType i = 0; i < 64; ++i) {
          word <<= 1;
          word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
        }
        _bv[level][cur_pos >> 6] = word;
      }
      if (size & 63ULL) {
        uint64_t word = 0ULL;
        for (SizeType i = 0; i < size - cur_pos; ++i) {
          word <<= 1;
          word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        _bv[level][size >> 6] = word;
      }
    }
  }

  auto get_bv() const {
    return _bv;
  }

private:
  std::vector<uint64_t*> _bv;
}; // class wt_ps

#endif // WT_PREFIX_SORTING

/******************************************************************************/