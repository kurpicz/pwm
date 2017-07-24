/*******************************************************************************
 * include/wt_ppc.hpp
 *
 * Copyright (C) 2016 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WT_PREFIX_COUNTING_PARALLEL
#define WT_PREFIX_COUNTING_PARALLEL

#include <cstring>
#include <omp.h>
#include <vector>

template <typename AlphabetType>
class wt_ppc {

public:
  static constexpr bool    is_parallel = true;
  static constexpr bool    is_tree     = true;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);

  wt_ppc() = default;

  wt_ppc(const std::vector<AlphabetType>& text, const uint64_t size,
    const uint64_t levels) : _bv(levels) {

    if(text.size() == 0) { return; }

    for (uint64_t level = 0; level < levels; ++level) {
      _bv[level] = new uint64_t[(size + 63ULL) >> 6];
      memset(_bv[level], 0, ((size + 63ULL) >> 6) * sizeof(uint64_t));
    }

    std::vector<uint64_t> hist;
    #pragma omp single
    {
      hist = std::vector<uint64_t>((1 << levels) * levels, 0);
    }

    #pragma omp barrier

    uint64_t* const hist_data_ptr = hist.data();

    #pragma omp parallel num_threads(levels) firstprivate(hist_data_ptr)
    {
      const uint64_t omp_rank = uint64_t(omp_get_thread_num());
      const uint64_t omp_size = uint64_t(omp_get_num_threads());

      const uint64_t global_max_char = (1 << levels);

      const AlphabetType* const text_ptr = text.data();

      uint64_t* const hist_ptr = hist_data_ptr + (global_max_char * omp_rank);

      #pragma omp for
      for (uint64_t cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < 64; ++i) {
          ++hist_ptr[text_ptr[cur_pos + i]];
          word <<= 1;
          word |= ((text_ptr[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        _bv[0][cur_pos >> 6] = word;
      }

      if ((size & 63ULL) && omp_rank == 0) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < (size & 63ULL); ++i) {
          ++hist[text[size - (size & 63ULL) + i]];
          word <<= 1;
          word |= ((text[size - (size & 63ULL) + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        _bv[0][size >> 6] = word;
      }
      #pragma omp for
      for (uint64_t i = 0; i < global_max_char; ++i) {
        for (uint64_t rank = 1; rank < omp_size; ++rank) {
          hist_data_ptr[i] += hist_data_ptr[(rank * global_max_char) + i];
        }
      }

      #pragma omp for
      for (uint64_t level = 1; level < levels; ++level) {
        const uint64_t local_max_char = (1 << level);
        const uint64_t requierd_characters = (1 << (levels - level));
        uint64_t* const hist_ttt = hist.data() + (level * global_max_char);
        for (uint64_t i = 0; i < local_max_char; ++i) {
          hist_ttt[i] = 0;
          for (uint64_t j = 0; j < requierd_characters; ++j) {
            hist_ttt[i] += hist_data_ptr[(i * requierd_characters) + j];
          }
        }
      }

      #pragma omp for
      for (uint64_t level = 1; level < levels; ++level) {

        const uint64_t local_max_char = (1 << level);
        const uint64_t prefix_shift = (levels - level);
        const uint64_t cur_bit_shift = prefix_shift - 1;

        uint64_t* const hist_ttt = hist.data() + (level * global_max_char);
        std::vector<uint64_t> borders(local_max_char);

        borders[0] = 0;
        for (uint64_t i = 1; i < local_max_char; ++i) {
          borders[i] = (borders[i - 1] + hist_ttt[i - 1]);
        }

        for (uint64_t i = 0; i < size; ++i) {
          const uint64_t pos = borders[text_ptr[i] >> prefix_shift]++;
          _bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
            << (63ULL - (pos & 63ULL)));
        }
      }
    }
  }

  auto get_bv_and_zeros() const {
    return std::make_pair(_bv, std::vector<uint64_t>());
  }

  auto get_bv() const {
    return _bv;
  }

private:
  std::vector<uint64_t*> _bv;
}; // class wt_ppc

#endif // WT_PREFIX_COUNTING_PARALLEL

/******************************************************************************/
