/*******************************************************************************
 * include/wm_cilk_ppc.hpp
 *
 * Copyright (C) 2016 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_CILK_PREFIX_COUNTING_PARALLEL
#define WM_CILK_PREFIX_COUNTING_PARALLEL

#include <cstring>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <vector>

#include "common.hpp"

template <typename TextType, typename SizeType>
class wm_cilk_ppc {

public:
  wm_cilk_ppc(const std::vector<TextType>& text, const SizeType size,
    const SizeType levels) : _bv(levels), _zeros(levels, 0) {
    
    for (SizeType level = 0; level < levels; ++level) {
      _bv[level] = new uint64_t[(size + 63ULL) >> 6];
        memset(_bv[level], 0, ((size + 63ULL) >> 6) * sizeof(uint64_t));  
    }

    std::vector<SizeType> hist;
    // #pragma single
    // {
      hist = std::vector<SizeType>((1 << levels) * levels, 0);
    // }

    // #pragma omp barrier

    SizeType* const hist_data_ptr = hist.data();

    // #pragma omp parallel num_threads(levels) firstprivate(hist_data_ptr)
    // {
      const auto omp_rank = __cilkrts_get_worker_number();
      const auto omp_size = __cilkrts_get_nworkers();

      const SizeType global_max_char = (1 << levels);

      const TextType* const text_ptr = text.data();

      SizeType* const hist_ptr = hist_data_ptr + (global_max_char * omp_rank);
      
      // While initializing the histogram, we also compute the fist level
      // #pragma omp for
      cilk_for (SizeType cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (SizeType i = 0; i < 64; ++i) {
          ++hist_ptr[text_ptr[cur_pos + i]];
          word <<= 1;
          word |= ((text_ptr[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        _bv[0][cur_pos >> 6] = word;
      }

      if ((size & 63ULL) && omp_rank == 0) {
        uint64_t word = 0ULL;
        for (SizeType i = 0; i < (size & 63ULL); ++i) {
          ++hist[text[size - (size & 63ULL) + i]];
          word <<= 1;
          word |= ((text[size - (size & 63ULL) + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (63 - (size & 63ULL));
        _bv[0][size >> 6] = word;
      }

      // Compute the historam with respect to the local slices of the text
      // #pragma omp for
      cilk_for (SizeType i = 0; i < global_max_char; ++i) {
        for (SizeType rank = 1; rank < omp_size; ++rank) {
          hist_data_ptr[i] += hist_data_ptr[(rank * global_max_char) + i];
        }
      }

      // The number of 0s at the last level is the number of "even" characters
      // #pragma omp single
      cilk_for (SizeType i = 0; i < global_max_char; i += 2) {
        _zeros[levels - 1] += hist_data_ptr[i];
      }
      // _zeros[levels - 1] = __sec_reduce_add (hist_data_ptr[0:global_max_char:2]);

      // Compute the histogram for each level of the WM.
      // #pragma omp for
      cilk_for (SizeType level = 1; level < levels; ++level) {
        const SizeType local_max_char = (1 << level);
        const SizeType requierd_characters = (1 << (levels - level));
        SizeType* const hist_ttt = hist.data() + (level * global_max_char);
        for (SizeType i = 0; i < local_max_char; ++i) {
          hist_ttt[i] = 0;
          for (SizeType j = 0; j < requierd_characters; ++j) {
            hist_ttt[i] += hist_data_ptr[(i * requierd_characters) + j];
          }
        }
      }

      // Now we compute the WM bottom-up, i.e., the last level first
      // #pragma omp for
      cilk_for (SizeType level = 1; level < levels; ++level) {

        const SizeType local_max_char = (1 << level);
        const SizeType prefix_shift = (levels - level);
        const SizeType cur_bit_shift = prefix_shift - 1;

        SizeType* const hist_ttt = hist.data() + (level * global_max_char);
        std::vector<SizeType> borders(local_max_char);
        std::vector<SizeType> bit_reverse = BitReverse<SizeType>(level);

        // Compute the starting positions of characters with respect to their
        // bit prefixes and the bit-reversal permutation
        borders[0] = 0;
        for (SizeType i = 1; i < local_max_char; ++i) {
          borders[bit_reverse[i]] = (borders[bit_reverse[i - 1]] +
             hist_ttt[bit_reverse[i - 1]]);
        }
        // The number of 0s is the position of the first 1 in the previous level
        _zeros[level - 1] = borders[1];

        // Now we insert the bits with respect to their bit prefixes
        for (SizeType i = 0; i < size; ++i) {
          const SizeType pos = borders[text_ptr[i] >> prefix_shift]++;
          _bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
            << (63ULL - (pos & 63ULL)));
        }
      }
    // }
  }

  auto get_bv_and_zeros() const {
    return std::make_pair(_bv, _zeros);
  }

private:
  std::vector<uint64_t*> _bv;
  std::vector<SizeType> _zeros;
}; // class wm_cilk_ppc

#endif // WM_CILK_PREFIX_COUNTING_PARALLEL

/******************************************************************************/