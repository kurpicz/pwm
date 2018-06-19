/*******************************************************************************
 * include/wx_ppc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <omp.h>

#include "util/common.hpp"
#include "util/ctx_all_levels.hpp"
#include "util/wavelet_structure.hpp"

template <typename AlphabetType, bool is_tree_, bool is_semi_external = false>
class wx_ppc {

public:
  static constexpr bool    is_parallel = true;
  static constexpr bool    is_tree     = is_tree_;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);
  static constexpr bool  is_huffman_shaped = false;

  using ctx_t = ctx_all_levels<is_tree, is_semi_external>;

  template <typename InputType>
  static wavelet_structure<is_semi_external> compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    if (size == 0) {
      return wavelet_structure<is_semi_external>();
    }

    const auto rho = rho_dispatch<is_tree>::create(levels);
    // TODO create new context with all max size hist-levels.
    ctx_t ctx(size, levels, rho);
    auto& bv = ctx.bv().raw_data();
    auto& zeros = ctx.zeros();

    std::vector<uint64_t> initial_hist((1ULL << levels) * levels, 0);

    #pragma omp parallel num_threads(levels)
    {
      const uint64_t omp_rank = uint64_t(omp_get_thread_num());
      const uint64_t omp_size = uint64_t(omp_get_num_threads());
      const uint64_t max_char = (1 << levels);

      auto* const initial_hist_ptr = 
        initial_hist.data() + (max_char * omp_rank);

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

      if ((size & 63ULL) && omp_rank == 0) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < (size & 63ULL); ++i) {
          ++initial_hist_ptr[text[size - (size & 63ULL) + i]];
          word <<= 1;
          word |= ((text[size - (size & 63ULL) + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        bv[0][size >> 6] = word;
      }

      #pragma omp barrier

      // Compute the historam with respect to the local slices of the text
      #pragma omp for
      for (uint64_t i = 0; i < max_char; ++i) {
        for (uint64_t rank = 0; rank < omp_size; ++rank) {
          ctx.hist(levels, i) += *(initial_hist.data() + (max_char * rank) + i);
        }
      }

      #pragma omp single
      {
        if (ctx_t::compute_zeros) {
          // The number of 0s at the last level is the number of "even" characters
          for (uint64_t i = 0; i < max_char; i += 2) {
            zeros[levels - 1] += ctx.hist(levels, i);
          }
        }
      }

      // Compute the histogram for each level of the wavelet structure
      #pragma omp for
      for (uint64_t level = 1; level < levels; ++level) {
        const uint64_t local_max_char = (1 << level);
        const uint64_t requierd_characters = (1 << (levels - level));
        for (uint64_t i = 0; i < local_max_char; ++i) {
          for (uint64_t j = 0; j < requierd_characters; ++j) {
            ctx.hist(level, i) +=
              ctx.hist(levels, (i * requierd_characters) + j);
          }
        }
      }

      // Now we compute the wavelet structure bottom-up, i.e., the last level first
      #pragma omp for
      for (uint64_t level = 1; level < levels; ++level) {

        const uint64_t local_max_char = (1 << level);
        const uint64_t prefix_shift = (levels - level);
        const uint64_t cur_bit_shift = prefix_shift - 1;

        // TODO: Add this local borders to the "new" context, too.
        std::vector<uint64_t> borders(local_max_char, 0);

        borders[0] = 0;
        for (uint64_t i = 1; i < local_max_char; ++i) {
          const auto prev_rho = ctx.rho(level, i - 1);
          borders[ctx.rho(level, i)] =
            borders[prev_rho] + ctx.hist(level, prev_rho);
        }
        // The number of 0s is the position of the first 1 in the previous level
        if (ctx_t::compute_zeros) {
          zeros[level - 1] = borders[1];
        }

        // Now we insert the bits with respect to their bit prefixes
        for (uint64_t i = 0; i < size; ++i) {
          const uint64_t pos = borders[text[i] >> prefix_shift]++;
          bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
            << (63ULL - (pos & 63ULL)));
        }
      }
    }
    if (ctx_t::compute_zeros) {
      return wavelet_structure<is_semi_external>(std::move(ctx.bv()), std::move(zeros));
    } else {
      return wavelet_structure<is_semi_external>(std::move(ctx.bv()));
    }
  }
}; // class wx_ppc

/******************************************************************************/
