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

#include <omp.h>
#include <vector>

#include <common.hpp>

template <typename AlphabetType, typename SizeType>
class wt_pps {

public:
  wt_pps(const std::vector<AlphabetType>& text, const SizeType size,
    const SizeType levels) : _bv(levels) {

    std::vector<SizeType*> borders;
    std::vector<SizeType*> hist;
    std::vector<AlphabetType> sorted_text(size);
    std::vector<SizeType> offsets(1 << levels, 0);

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

      for (SizeType level = levels - 1; level > 0; --level) {
        const SizeType prefix_shift = (levels - level);
        const SizeType cur_bit_shift = prefix_shift - 1;

        #pragma omp for
        for (int32_t rank = 0; rank < omp_size; ++rank) {
          for (SizeType i = 0;
               i < global_max_char; i += (1ULL << prefix_shift)) {
            hist[rank][i] += hist[rank][i + (1ULL << cur_bit_shift)];
          }
        }

        #pragma omp for
        for (SizeType i = 0; i < global_max_char; i += (1ULL << prefix_shift)) {
          borders[0][i] = 0;
          for (int32_t rank = 1; rank < omp_size; ++rank) {
            borders[rank][i] = borders[rank - 1][i] + hist[rank - 1][i];
          }
        }

        #pragma omp single
        {
          for (SizeType i = 1; i < (1ULL << level); ++i) {
            offsets[i << prefix_shift] = 
              offsets[(i - 1) << prefix_shift] + 
              borders[omp_size - 1][(i - 1) << prefix_shift] +
              hist[omp_size - 1][(i - 1) << prefix_shift];
          }
        }

        #pragma omp for
        for (SizeType i = 0; i <= size - 64; i += 64) {
          for (SizeType j = 0; j < 64; ++j) {
            const AlphabetType considerd_char =
              (text[i + j] >> prefix_shift) << (prefix_shift);
            const SizeType bucket_border = borders[omp_rank][considerd_char]++;
            const SizeType offset = offsets[considerd_char];
            sorted_text[bucket_border + offset] = text[i + j];
          }
        }
        if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
          for (SizeType i = size - (size & 63ULL); i < size; ++i) {
            const AlphabetType considerd_char =
              (text[i] >> prefix_shift) << (prefix_shift);
            const SizeType bucket_border = borders[omp_rank][considerd_char]++;
            const SizeType offset = offsets[considerd_char];
            sorted_text[offset + bucket_border] = text[i];
          }
        }

        #pragma omp for
        for (SizeType cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
          uint64_t word = 0ULL;
          for (SizeType i = 0; i < 64; ++i) {
            word <<= 1;
            word |= ((sorted_text[cur_pos + i] >> cur_bit_shift) & 1ULL);
          }
          _bv[level][cur_pos >> 6] = word;
        }
        if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
          uint64_t word = 0ULL;
          for (SizeType i = 0; i < (size & 63ULL); ++i) {
            word <<= 1;
            word |= ((sorted_text[size - (size & 63ULL) + i]
              >> cur_bit_shift) & 1ULL);
          }
          word <<= (63 - (size & 63ULL));
          _bv[level][size >> 6] = word;
        }
      }
    }
  }

  auto get_bv() const {
    return _bv;
  }

private:
  std::vector<uint64_t*> _bv;
}; // class wt_pps

#endif // WT_PREFIX_SORTING_PARALLEL

/******************************************************************************/