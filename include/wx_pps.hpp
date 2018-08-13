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
#include "util/ctx_sliced_single_level.hpp"
#include "util/wavelet_structure.hpp"
#include "util/memory_types.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_pps {

public:
  static constexpr bool    is_parallel = true;
  static constexpr bool    is_tree     = is_tree_;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);
  static constexpr bool  is_huffman_shaped = false;
  static constexpr memory_mode mem_mode = memory_mode::internal;

  template <typename InputType>
  static wavelet_structure compute(const InputType& text, const uint64_t size,
    const uint64_t levels) {

    using ctx_t = ctx_sliced_single_level<is_tree>;

    if (size == 0) {
      return wavelet_structure();
    }

    ctx_t ctx;
    std::vector<AlphabetType> sorted_text(size);
    std::vector<uint64_t> offsets(1 << levels, 0);

    #pragma omp parallel
    {
      const auto omp_rank = omp_get_thread_num();
      const auto omp_size = omp_get_num_threads();
      const uint64_t max_char = (1 << levels);

      #pragma omp single
      ctx = ctx_t(size, levels, omp_size);

      auto& bv = ctx.bv().raw_data();
      auto& zeros = ctx.zeros();

      // While initializing the histogram, we also compute the fist level
      #pragma omp for
      for (uint64_t cur_pos = 0; cur_pos <= size - 64; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < 64; ++i) {
          ++ctx.hist(omp_rank, text[cur_pos + i]);
          word <<= 1;
          word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        bv[0][cur_pos >> 6] = word;
      }
      if ((size & 63ULL) && ((omp_rank + 1) == omp_size)) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < (size & 63ULL); ++i) {
          ++ctx.hist(omp_rank, text[size - (size & 63ULL) + i]);
          word <<= 1;
          word |= ((text[size - (size & 63ULL) + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        bv[0][size >> 6] = word;
      }

      // The number of 0's at the last level is the number of "even" characters
      #pragma omp single
      for (uint64_t i = 0; i < max_char; i += 2) {
        for (int32_t rank = 0; rank < omp_size; ++rank) {
          zeros[levels - 1] += ctx.hist(rank, i);
        }
      }

      // Now we compute the wavelet structure bottom-up, i.e., the last level first
      for (uint64_t level = levels - 1; level > 0; --level) {
        const uint64_t prefix_shift = (levels - level);
        const uint64_t cur_bit_shift = prefix_shift - 1;

        // Compute the histogram and the border for each bit prefix and
        // processor, i.e., for one fixed bit prefix we compute the prefix sum
        // over the number of occurrences at each processor
        #pragma omp for
        for (uint64_t i = 0; i < max_char; i += (1ULL << prefix_shift)) {
          ctx.borders(0, i) = 0;
          ctx.hist(0, i) += ctx.hist(0, i + (1ULL << cur_bit_shift));
          for (int32_t rank = 1; rank < omp_size; ++rank) {
            ctx.hist(rank, i) += ctx.hist(rank, i + (1ULL << cur_bit_shift));
            ctx.borders(rank, i) = 
              ctx.borders(rank - 1, i) + ctx.hist(rank - 1, i);
          }
        }

        // Now we compute the offset for each bit prefix, i.e., the number of
        // lexicographically smaller characters
        #pragma omp single
        {
          for (uint64_t i = 1; i < (1ULL << level); ++i) {
            const auto prev_rho = ctx.rho(level, i - 1);
            offsets[ctx.rho(level, i) << prefix_shift] =
              offsets[ctx.rho(level, i - 1) << prefix_shift] +
              ctx.borders(omp_size - 1, prev_rho << prefix_shift) +
              ctx.hist(omp_size - 1, prev_rho << prefix_shift);
            if (ctx_t::compute_rho)  {
              ctx.set_rho(level, i - 1, prev_rho >> 1);
            }
          }
          // The number of 0s is the position of the first 1 at the first
          // processor
          if (ctx_t::compute_zeros) {
            zeros[level - 1] = offsets[1ULL << prefix_shift];
          }
        }
        // We add the offset to the borders (for performance)
        #pragma omp for
        for (int32_t rank = 0; rank < omp_size; ++rank) {
          for (uint64_t i = 0; i < max_char; i += (1ULL << prefix_shift)) {
            ctx.borders(rank, i) += offsets[i];
          }
        }

        // We align the borders (in memory) to increase performance by reducing
        // the number of cache misses
        std::vector<uint64_t> borders_aligned(1ULL << level, 0);
        for (uint64_t i = 0; i < max_char; i += (1ULL << prefix_shift)) {
          borders_aligned[i >> prefix_shift] = ctx.borders(omp_rank, i);
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

        #pragma omp barrier

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
    if (ctx_t::compute_zeros) {
      return wavelet_structure(std::move(ctx.bv()), std::move(ctx.zeros()));
    } else {
      return wavelet_structure(std::move(ctx.bv()));
    }
  }
};

/******************************************************************************/
