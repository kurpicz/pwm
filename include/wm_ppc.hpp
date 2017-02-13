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

#include <vector>
#include <omp.h>

#include "common.hpp"

template <typename TextType, typename SizeType>
class wm_ppc {

public:
  wm_ppc(const std::vector<TextType>& text, const SizeType size,
    const SizeType levels) : _bv(levels), _zeros(levels, 0) {
    
    for (SizeType level = 0; level < levels; ++level) {
      _bv[level] = new uint64_t[(size + 63ULL) >> 6];
        memset(_bv[level], 0, ((size + 63ULL) >> 3));  
    }

    std::vector<SizeType> hist;
    #pragma single
    {
      hist = std::vector<SizeType>((1 << levels) * levels, 0);
    }

    #pragma omp barrier

    SizeType* const hist_data_ptr = hist.data();

    #pragma omp parallel num_threads(levels) firstprivate(hist_data_ptr)
    {
      const auto omp_rank = omp_get_thread_num();
      const auto omp_size = omp_get_num_threads();

      const SizeType global_max_char = (1 << levels);

      const TextType* const text_ptr = text.data();

      SizeType* const hist_ptr = hist_data_ptr + (global_max_char * omp_rank);
      
      #pragma omp for
      for (SizeType cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
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
      #pragma omp for
      for (SizeType i = 0; i < global_max_char; ++i) {
        for (SizeType rank = 1; rank < omp_size; ++rank) {
          hist_data_ptr[i] += hist_data_ptr[(rank * global_max_char) + i];
        }
      }

      #pragma omp single
      for (SizeType i = 0; i < global_max_char; i += 2) {
      	_zeros[levels - 1] += hist_data_ptr[i];
      }

      #pragma omp for
      for (SizeType level = 1; level < levels; ++level) {
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

      #pragma omp for
      for (SizeType level = 1; level < levels; ++level) {

        const SizeType local_max_char = (1 << level);
        const SizeType prefix_shift = (levels - level);
        const SizeType cur_bit_shift = prefix_shift - 1;

        SizeType* const hist_ttt = hist.data() + (level * global_max_char);
        std::vector<SizeType> borders(local_max_char);
        std::vector<SizeType> bit_reverse = BitReverse<SizeType>(level);

        borders[0] = 0;
        for (SizeType i = 1; i < local_max_char; ++i) {
          borders[bit_reverse[i]] = (borders[bit_reverse[i - 1]] +
             hist_ttt[bit_reverse[i - 1]]);
        }
        _zeros[level - 1] = borders[1];

        for (SizeType i = 0; i < size; ++i) {
          const SizeType pos = borders[text_ptr[i] >> prefix_shift]++;
          _bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
            << (63ULL - (pos & 63ULL)));
        }
      }
    }
  }

  ~wm_ppc() {
    for (uint64_t i = 0; i < _bv.size(); ++i) {
      delete [] _bv[i];
    }
  }

  auto get_bv_and_zeros() const {
    return std::make_pair(_bv, _zeros);
  }

private:
  std::vector<uint64_t*> _bv;
  std::vector<SizeType> _zeros;
}; // class wm_ppc

#endif // WM_PREFIX_COUNTING_PARALLEL

/******************************************************************************/