/*******************************************************************************
 * include/wm_prefix_counting.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_PREFIX_COUNTING
#define WM_PREFIX_COUNTING

#include <cstring>
#include <vector>

#include "common.hpp"

template <typename AlphabetType, typename SizeType>
class wm_pc {

public:
  wm_pc(const std::vector<AlphabetType>& text, const SizeType size,
    const SizeType levels) : _bv(levels), _zeros(levels, 0) {

    SizeType cur_max_char = (1 << levels);
    std::vector<SizeType> bit_reverse = BitReverse<SizeType>(levels - 1);
    std::vector<SizeType> s_pos(cur_max_char, 0);
    std::vector<SizeType> hist(cur_max_char, 0);
    std::vector<SizeType> borders(cur_max_char, 0);

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
      _bv[0][cur_pos >> 6] = word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (SizeType i = 0; i < size - cur_pos; ++i) {
        ++hist[text[cur_pos + i]];
        word <<= 1;
        word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      }
      word <<= (64 - (size & 63ULL));
      _bv[0][size >> 6] = word;
    }

    for (SizeType i = 0; i < cur_max_char; i += 2) {
      _zeros[levels - 1] += hist[i];
    }

    for (SizeType level = levels - 1; level > 0; --level) {
      const SizeType prefix_shift = (levels - level);
      const SizeType cur_bit_shift = prefix_shift - 1;

      _bv[level] = new uint64_t[(size + 63ULL) >> 6];
      memset(_bv[level], 0, ((size + 63ULL) >> 3)); // memset is ok (all to 0)

      cur_max_char >>= 1;
      for (SizeType i = 0; i < cur_max_char; ++i) {
        hist[i] = hist[i << 1] + hist[(i << 1) + 1];
      }

      borders[0] = 0;
      for (SizeType i = 1; i < cur_max_char; ++i) {
        borders[bit_reverse[i]] = borders[bit_reverse[i - 1]] +
          hist[bit_reverse[i - 1]];
        bit_reverse[i - 1] >>= 1;
      }
      _zeros[level - 1] = borders[1];

      for (SizeType i = 0; i < size; ++i) {
        const SizeType pos = borders[text[i] >> prefix_shift]++;
        _bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
          << (63ULL - (pos & 63ULL)));
      }
    }
  }

  auto get_bv_and_zeros() const {
    return std::make_pair(_bv, _zeros);
  }

private:
  std::vector<uint64_t*> _bv;
  std::vector<SizeType> _zeros;
}; // class wm_pc

#endif // WM_PREFIX_COUNTING

/******************************************************************************/