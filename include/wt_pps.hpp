/*******************************************************************************
 * include/wt_pps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WT_PREFIX_SORTING_PARALLEL
#define WT_PREFIX_SORTING_PARALLEL

#include <cstring>
#include <omp.h>
#include <vector>

template <typename AlphabetType>
class wt_pps {

public:
  static constexpr bool    is_parallel = true;
  static constexpr bool    is_tree     = true;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);

    static wavelet_structure compute(AlphabetType const* const text,
                                     const uint64_t size,
                                     const uint64_t levels)
    {

    if(size == 0) { return wavelet_structure(); }

    auto _bv = Bvs(size, levels);
    auto& bv = _bv.vec();

    std::vector<uint64_t*> borders;
    std::vector<uint64_t*> hist;
    std::vector<AlphabetType> sorted_text(size);
    std::vector<uint64_t> offsets(1 << levels, 0);

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
      const uint64_t omp_rank = uint64_t(omp_get_thread_num());
      const uint64_t omp_size = uint64_t(omp_get_num_threads());
      const uint64_t global_max_char = (1 << levels);

      #pragma omp for
      for (uint64_t cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < 64; ++i) {
          ++hist[omp_rank][text[cur_pos + i]];
          word <<= 1;
          word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        bv[0][cur_pos >> 6] = word;
      }
      if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < (size & 63ULL); ++i) {
          ++hist[omp_rank][text[size - (size & 63ULL) + i]];
          word <<= 1;
          word |= ((text[size - (size & 63ULL) + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        bv[0][size >> 6] = word;
      }

      for (uint64_t level = levels - 1; level > 0; --level) {
        const uint64_t prefix_shift = (levels - level);
        const uint64_t cur_bit_shift = prefix_shift - 1;

        #pragma omp for
        for (uint64_t i = 0; i < global_max_char; i += (1ULL << prefix_shift)) {
          borders[0][i] = 0;
          hist[0][i] += hist[0][i + (1ULL << cur_bit_shift)];
          for (uint64_t rank = 1; rank < omp_size; ++rank) {
            hist[rank][i] += hist[rank][i + (1ULL << cur_bit_shift)];
            borders[rank][i] = borders[rank - 1][i] + hist[rank - 1][i];
          }
        }

        #pragma omp single
        {
          for (uint64_t i = 1; i < (1ULL << level); ++i) {
            offsets[i << prefix_shift] =
              offsets[(i - 1) << prefix_shift] +
              borders[omp_size - 1][(i - 1) << prefix_shift] +
              hist[omp_size - 1][(i - 1) << prefix_shift];
          }
        }

        #pragma omp for
        for (uint64_t rank = 0; rank < omp_size; ++rank) {
          for (uint64_t i = 0; i < global_max_char; i += (1ULL << prefix_shift)) {
            borders[rank][i] += offsets[i];
          }
        }

        std::vector<uint64_t> borders_aligned(1ULL << level, 0);
        #pragma omp simd
        for (uint64_t i = 0; i < global_max_char; i += (1ULL << prefix_shift)) {
          borders_aligned[i >> prefix_shift] = borders[omp_rank][i];
        }

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

        #pragma omp for
        for (uint64_t cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
          uint64_t word = 0ULL;
          for (uint64_t i = 0; i < 64; ++i) {
            word <<= 1;
            word |= (sorted_text[cur_pos + i] & 1ULL);
          }
          bv[level][cur_pos >> 6] = word;
        }
        if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
          uint64_t word = 0ULL;
          for (uint64_t i = 0; i < (size & 63ULL); ++i) {
            word <<= 1;
            word |= (sorted_text[size - (size & 63ULL) + i] & 1ULL);
          }
          word <<= (64 - (size & 63ULL));
          bv[level][size >> 6] = word;
        }
      }
    }
    return wavelet_structure(std::move(_bv));
  }
}; // class wt_pps

#endif // WT_PREFIX_SORTING_PARALLEL

/******************************************************************************/
