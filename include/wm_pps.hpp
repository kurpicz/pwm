/*******************************************************************************
 * include/wm_pps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_PREFIX_SORTING_PARALLEL
#define WM_PREFIX_SORTING_PARALLEL

#include <cstring>
#include <omp.h>
#include <vector>

#include "util/common.hpp"

template <typename AlphabetType>
class wm_pps {

public:
  static constexpr bool    is_parallel = true;
  static constexpr bool    is_tree     = false;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);

  wm_pps() = default;

  wm_pps(const std::vector<AlphabetType>& text, const uint64_t size,
    const uint64_t levels) : _bv(levels), _zeros(levels, 0) {

    if(text.size() == 0) { return; }

    std::vector<uint64_t*> borders;
    std::vector<uint64_t*> hist;
    std::vector<AlphabetType> sorted_text(size);
    std::vector<uint64_t> offsets(1 << levels, 0);
    std::vector<uint64_t> bit_reverse = BitReverse(levels - 1);

    for (uint64_t level = 0; level < levels; ++level) {
      _bv[level] = new uint64_t[(size + 63ULL) >> 6];
      memset(_bv[level], 0, ((size + 63ULL) >> 6) * sizeof(uint64_t));
    }

    int32_t num_threads;
    #pragma omp parallel
    {
      num_threads = omp_get_num_threads();
    }

    #pragma omp single
    {
      hist.reserve(num_threads);
      borders.reserve(num_threads);
      for (int32_t rank = 0; rank < num_threads; ++rank) {
        hist[rank] = new uint64_t[1 << levels];
        memset(hist[rank], 0, (1 << levels) * sizeof(uint64_t));
        borders[rank] = new uint64_t[1 << levels];
        memset(borders[rank], 0, (1 << levels) * sizeof(uint64_t));
      }
    }

    #pragma omp parallel
    {
      const auto omp_rank = omp_get_thread_num();
      const auto omp_size = omp_get_num_threads();

      const uint64_t global_max_char = (1 << levels);
      uint64_t cur_max_char = global_max_char;

      // While initializing the histogram, we also compute the fist level
      #pragma omp for
      for (uint64_t cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < 64; ++i) {
          ++hist[omp_rank][text[cur_pos + i]];
          word <<= 1;
          word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        _bv[0][cur_pos >> 6] = word;
      }
      if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < (size & 63ULL); ++i) {
          ++hist[omp_rank][text[size - (size & 63ULL) + i]];
          word <<= 1;
          word |= ((text[size - (size & 63ULL) + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        _bv[0][size >> 6] = word;
      }

      // The number of 0s at the last level is the number of "even" characters
      #pragma omp single
      for (uint64_t i = 0; i < cur_max_char; i += 2) {
        for (int32_t rank = 0; rank < omp_size; ++rank) {
          _zeros[levels - 1] += hist[rank][i];
        }
      }

      // Now we compute the WM bottom-up, i.e., the last level first
      for (uint64_t level = levels - 1; level > 0; --level) {
        const uint64_t prefix_shift = (levels - level);
        const uint64_t cur_bit_shift = prefix_shift - 1;

        // Compute the histogram and the border for each bit prefix and
        // processor, i.e., for one fixed bit prefix we compute the prefix sum
        // over the number of occurrences at each processor
        #pragma omp for
        for (uint64_t i = 0; i < global_max_char; i += (1ULL << prefix_shift)) {
          borders[0][i] = 0;
          hist[0][i] += hist[0][i + (1ULL << cur_bit_shift)];
          for (int32_t rank = 1; rank < omp_size; ++rank) {
            hist[rank][i] += hist[rank][i + (1ULL << cur_bit_shift)];
            borders[rank][i] = borders[rank - 1][i] + hist[rank - 1][i];
          }
        }

        // Now we compute the offset for each bit prefix, i.e., the number of
        // lexicographically smaller characters
        #pragma omp single
        {
          for (uint64_t i = 1; i < (1ULL << level); ++i) {
            offsets[bit_reverse[i] << prefix_shift] =
              offsets[bit_reverse[i - 1] << prefix_shift] +
              borders[omp_size - 1][bit_reverse[i - 1] << prefix_shift] +
              hist[omp_size - 1][bit_reverse[i - 1] << prefix_shift];
            bit_reverse[i - 1] >>= 1;
          }
          // The number of 0s is the position of the first 1 at the first
          // processor
          _zeros[level - 1] = offsets[1ULL << prefix_shift];
        }

        // We add the offset to the borders (for performance)
        #pragma omp for
        for (int32_t rank = 0; rank < omp_size; ++rank) {
          for (uint64_t i = 0; i < global_max_char; i += (1ULL << prefix_shift)) {
            borders[rank][i] += offsets[i];
          }
        }

        // We align the borders (in memory) to increase performance by reducing
        // the number of cache misses
        std::vector<uint64_t> borders_aligned(1ULL << level, 0);
        #pragma omp simd
        for (uint64_t i = 0; i < global_max_char; i += (1ULL << prefix_shift)) {
          borders_aligned[i >> prefix_shift] = borders[omp_rank][i];
        }

        // Sort the text using the computed (and aligned) borders
        #pragma omp for
        for (uint64_t i = 0; i <= size - 64; i += 64) {
          for (uint64_t j = 0; j < 64; ++j) {
            const AlphabetType considerd_char = (text[i + j] >> cur_bit_shift);
            sorted_text[borders_aligned[considerd_char >> 1]++] = considerd_char;
          }
        }
        if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
          for (uint64_t i = size - (size & 63ULL); i < size; ++i) {
            const AlphabetType considerd_char = (text[i] >> cur_bit_shift);
            sorted_text[borders_aligned[considerd_char >> 1]++] = considerd_char;
          }
        }

        // Since we have sorted the text, we can simply scan it from left to
        // right and for the character at position $i$ we set the $i$-th bit in
        // the bit vector accordingly
        #pragma omp for
        for (uint64_t cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
          uint64_t word = 0ULL;
          for (uint64_t i = 0; i < 64; ++i) {
            word <<= 1;
            word |= (sorted_text[cur_pos + i] & 1ULL);
          }
          _bv[level][cur_pos >> 6] = word;
        }

        if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
          uint64_t word = 0ULL;
          for (uint64_t i = 0; i < (size & 63ULL); ++i) {
            word <<= 1;
            word |= (sorted_text[size - (size & 63ULL) + i] & 1ULL);
          }
          word <<= (64 - (size & 63ULL));
          _bv[level][size >> 6] = word;
        }
      }
    }
  }

  auto get_bv_and_zeros() const {
    return std::make_pair(_bv, _zeros);
  }

private:
  std::vector<uint64_t*> _bv;
  std::vector<uint64_t> _zeros;
}; // class wm_pps

#endif // WM_PREFIX_SORTING_PARALLEL

/******************************************************************************/
