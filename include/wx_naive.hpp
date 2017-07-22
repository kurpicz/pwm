/*******************************************************************************
 * include/wx_naive.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>

template <typename AlphabetType, bool is_matrix>
class wx_naive;

template <typename AlphabetType>
class wx_naive<AlphabetType, false> {

public:
  wx_naive(const std::vector<AlphabetType>& text, const uint64_t size,
    const uint64_t levels) : _bv(levels) {

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

template <typename AlphabetType>
class wx_naive<AlphabetType, true> {

public:
  wx_naive(const std::vector<AlphabetType>& text, const uint64_t size,
    const uint64_t levels) : _bv(levels), _zeros(levels, 0) {

    if(text.size() == 0) { return; }

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
  std::vector<uint64_t> _zeros;
}; // class wx_naive<MATRIX>

/******************************************************************************/
