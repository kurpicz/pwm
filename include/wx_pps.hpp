/*******************************************************************************
 * include/wx_pps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <omp.h>
#include "util/common.hpp"
#include "util/context.hpp"
#include "util/wavelet_structure.hpp"

template <typename AlphabetType, bool is_matrix>
class wx_pps {
  using ctx_t =
    KeepLevel<is_matrix, typename rho_dispatch<is_matrix>::type>;

public:
  static constexpr bool    is_parallel = true;
  static constexpr bool    is_tree     = !is_matrix;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);

  static wavelet_structure compute(AlphabetType const* const text,
    const uint64_t size, const uint64_t levels) {

    if (size == 0) {
      return wavelet_structure();
    }

    const auto rho = rho_dispatch<is_matrix>::create(levels);
    // TODO create new context with all max size hist-levels.
    ctx_t ctx(size, levels, rho);
    auto& bv = ctx.bv().vec();
    auto& zeros = ctx.zeros();

    std::vector<uint64_t> initial_hist;

    std::vector<AlphabetType> sorted_text(size);    
    #pragma omp parallel
    {
      const auto omp_rank = omp_get_thread_num();
      const auto omp_size = omp_get_num_threads();
      const uint64_t max_char = (1 << levels);

      #pragma omp single
      initial_hist = std::vector<uint64_t>(max_char * omp_size, 0);

      auto* const initial_hist_ptr = 
        initial_hist.data() + (max_char * omp_rank);

      // While initializing the histogram, we also compute the fist level
      #pragma omp for
      for (uint64_t cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < 64; ++i) {
          ++initial_hist_ptr[text[cur_pos + i]];
          word <<= 1;
          word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        bv[0][cur_pos >> 6] = word;
      }

      if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < (size & 63ULL); ++i) {
          ++initial_hist_ptr[text[size - (size & 63ULL) + i]];
          word <<= 1;
          word |= ((text[size - (size & 63ULL) + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        bv[0][size >> 6] = word;
      }

    }
    if (ctx_t::compute_zeros) {
      return wavelet_structure(std::move(ctx.bv()), std::move(zeros));
    } else {
      return wavelet_structure(std::move(ctx.bv()));
    }
  }
};

/******************************************************************************/
