/*******************************************************************************
 * include/wm_ppc.hpp
 *
 * Copyright (C) 2016 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_PREFIX_COUNTING_PARALLEL
#define WM_PREFIX_COUNTING_PARALLEL

#include <cstring>
#include <omp.h>
#include <vector>

#include "util/common.hpp"

template <typename AlphabetType>
class wm_ppc {

public:
  static constexpr bool    is_parallel = true;
  static constexpr bool    is_tree     = false;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);

    static wavelet_structure compute(AlphabetType const* const text,
                                     const uint64_t size,
                                     const uint64_t levels)
    {

    if(size == 0) { return wavelet_structure(); }

    auto _zeros = std::vector<size_t>(levels, 0);
    auto _bv = Bvs(size, levels);
    auto& bv = _bv.vec();

    std::vector<uint64_t> hist;
    #pragma omp single
    {
      hist = std::vector<uint64_t>((1ULL << levels) * levels, 0);
    }

    #pragma omp barrier

    uint64_t* const hist_data_ptr = hist.data();

    #pragma omp parallel num_threads(levels) firstprivate(hist_data_ptr)
    {
      const uint64_t omp_rank = uint64_t(omp_get_thread_num());
      const uint64_t omp_size = uint64_t(omp_get_num_threads());

      const uint64_t global_max_char = (1ULL << levels);

      const AlphabetType* const text_ptr = text;

      uint64_t* const hist_ptr = hist_data_ptr + (global_max_char * omp_rank);

      // While initializing the histogram, we also compute the fist level
      #pragma omp for
      for (uint64_t cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < 64; ++i) {
          ++hist_ptr[text_ptr[cur_pos + i]];
          word <<= 1;
          word |= ((text_ptr[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        bv[0][cur_pos >> 6] = word;
      }

      if ((size & 63ULL) && omp_rank == 0) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < (size & 63ULL); ++i) {
          ++hist[text[size - (size & 63ULL) + i]];
          word <<= 1;
          word |= ((text[size - (size & 63ULL) + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        bv[0][size >> 6] = word;
      }

      // Compute the historam with respect to the local slices of the text
      #pragma omp for
      for (uint64_t i = 0; i < global_max_char; ++i) {
        for (uint64_t rank = 1; rank < omp_size; ++rank) {
          hist_data_ptr[i] += hist_data_ptr[(rank * global_max_char) + i];
        }
      }

      // The number of 0s at the last level is the number of "even" characters
      #pragma omp single
      for (uint64_t i = 0; i < global_max_char; i += 2) {
      	_zeros[levels - 1] += hist_data_ptr[i];
      }

      // Compute the histogram for each level of the WM.
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

      // Now we compute the WM bottom-up, i.e., the last level first
      #pragma omp for
      for (uint64_t level = 1; level < levels; ++level) {

        const uint64_t local_max_char = (1 << level);
        const uint64_t prefix_shift = (levels - level);
        const uint64_t cur_bit_shift = prefix_shift - 1;

        uint64_t* const hist_ttt = hist.data() + (level * global_max_char);
        std::vector<uint64_t> borders(local_max_char, 0);
        std::vector<uint64_t> bit_reverse = BitReverse(level);

        // Compute the starting positions of characters with respect to their
        // bit prefixes and the bit-reversal permutation
        borders[0] = 0;
        for (uint64_t i = 1; i < local_max_char; ++i) {
          borders[bit_reverse[i]] = (borders[bit_reverse[i - 1]] +
             hist_ttt[bit_reverse[i - 1]]);
        }
        // The number of 0s is the position of the first 1 in the previous level
        _zeros[level - 1] = borders[1];

        // Now we insert the bits with respect to their bit prefixes
        for (uint64_t i = 0; i < size; ++i) {
          const uint64_t pos = borders[text_ptr[i] >> prefix_shift]++;
          bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
            << (63ULL - (pos & 63ULL)));
        }
      }
    }

    return wavelet_structure(std::move(_bv), std::move(_zeros));
  }
}; // class wm_ppc

#endif // WM_PREFIX_COUNTING_PARALLEL

/******************************************************************************/
