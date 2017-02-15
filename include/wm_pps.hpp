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

#include <chrono>
#include <omp.h>
#include <vector>

#include <common.hpp>

template <typename AlphabetType, typename SizeType>
class wm_pps {

public:
  wm_pps(const std::vector<AlphabetType>& text, const SizeType size,
    const SizeType levels) : _bv(levels), _zeros(levels, 0) {

    std::vector<SizeType*> borders;
    std::vector<SizeType*> hist;
    std::vector<AlphabetType> sorted_text(size);
    std::vector<SizeType> offsets(1 << levels, 0);
    std::vector<SizeType> bit_reverse = BitReverse<SizeType>(levels - 1);

    for (SizeType level = 0; level < levels; ++level) {
      _bv[level] = new uint64_t[(size + 63ULL) >> 6];
      memset(_bv[level], 0, ((size + 63ULL) >> 3));
    }

    int32_t num_threads;
    #pragma omp parallel
    {
      num_threads = omp_get_num_threads();
    }

    #pragma single
    {
      hist.reserve(num_threads);
      borders.reserve(num_threads);
      for (int32_t rank = 0; rank < num_threads; ++rank) {
        hist[rank] = new SizeType[1 << levels];
        memset(hist[rank], 0, (1 << levels) * sizeof(SizeType));
        borders[rank] = new SizeType[1 << levels];
        memset(borders[rank], 0, (1 << levels) * sizeof(SizeType));
      }
    }

    #pragma omp parallel
    {
      const auto omp_rank = omp_get_thread_num();
      const auto omp_size = omp_get_num_threads();

      const SizeType global_max_char = (1 << levels);
      SizeType cur_max_char = global_max_char;

      #pragma omp for
      for (SizeType cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (SizeType i = 0; i < 64; ++i) {
          ++hist[omp_rank][text[cur_pos + i]];
          word <<= 1;
          word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        _bv[0][cur_pos >> 6] = word;
      }
      if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
        uint64_t word = 0ULL;
        for (SizeType i = 0; i < (size & 63ULL); ++i) {
          ++hist[omp_rank][text[size - (size & 63ULL) + i]];
          word <<= 1;
          word |= ((text[size - (size & 63ULL) + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (63 - (size & 63ULL));
        _bv[0][size >> 6] = word;
      }

      #pragma omp single
      for (SizeType i = 0; i < cur_max_char; i += 2) {
        for (int32_t rank = 0; rank < omp_size; ++rank) {
          _zeros[levels - 1] += hist[rank][i];
        }
      }

      for (SizeType level = levels - 1; level > 0; --level) {
        const SizeType prefix_shift = (levels - level);
        const SizeType cur_bit_shift = prefix_shift - 1;

        #pragma omp for
        for (SizeType i = 0; i < global_max_char; i += (1ULL << prefix_shift)) {
          borders[0][i] = 0;
          hist[0][i] += hist[0][i + (1ULL << cur_bit_shift)];
          for (int32_t rank = 1; rank < omp_size; ++rank) {
            hist[rank][i] += hist[rank][i + (1ULL << cur_bit_shift)];
            borders[rank][i] = borders[rank - 1][i] + hist[rank - 1][i];
          }
        }

        #pragma omp single
        {
          for (SizeType i = 1; i < (1ULL << level); ++i) {
            offsets[bit_reverse[i] << prefix_shift] = 
              offsets[bit_reverse[i - 1] << prefix_shift] + 
              borders[omp_size - 1][bit_reverse[i - 1] << prefix_shift] +
              hist[omp_size - 1][bit_reverse[i - 1] << prefix_shift];
            bit_reverse[i - 1] >>= 1;
          }
          _zeros[level - 1] = offsets[1ULL << prefix_shift];
        }

        #pragma omp for
        for (int32_t rank = 0; rank < omp_size; ++rank) {
          for (SizeType i = 0; i < global_max_char; i += (1ULL << prefix_shift)) {
            borders[rank][i] += offsets[i];
          }
        }

        std::vector<SizeType> borders_aligned(1ULL << level, 0);
        #pragma omp simd
        for (SizeType i = 0; i < global_max_char; i += (1ULL << prefix_shift)) {
          borders_aligned[i >> prefix_shift] = borders[omp_rank][i];
        }

        #pragma omp for
        for (SizeType i = 0; i <= size - 64; i += 64) {
          for (SizeType j = 0; j < 64; ++j) {
            const AlphabetType considerd_char = (text[i + j] >> cur_bit_shift);
            sorted_text[borders_aligned[considerd_char >> 1]++] = considerd_char;
          }
        }
        if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
          for (SizeType i = size - (size & 63ULL); i < size; ++i) {
            const AlphabetType considerd_char = (text[i] >> cur_bit_shift);
            sorted_text[borders_aligned[considerd_char >> 1]++] = considerd_char;
          }
        }

        #pragma omp for
        for (SizeType cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
          uint64_t word = 0ULL;
          for (SizeType i = 0; i < 64; ++i) {
            word <<= 1;
            word |= (sorted_text[cur_pos + i] & 1ULL);
          }
          _bv[level][cur_pos >> 6] = word;
        }
        if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
          uint64_t word = 0ULL;
          for (SizeType i = 0; i < (size & 63ULL); ++i) {
            word <<= 1;
            word |= (sorted_text[size - (size & 63ULL) + i] & 1ULL);
          }
          word <<= (63 - (size & 63ULL));
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
  std::vector<SizeType> _zeros;
}; // class wm_pps

#endif // WM_PREFIX_SORTING_PARALLEL

/******************************************************************************/