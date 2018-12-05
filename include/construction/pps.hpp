/*******************************************************************************
 * include/construction/pps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <omp.h>

#include "construction/building_blocks.hpp"
#include "construction/wavelet_structure.hpp"
#include "util/common.hpp"

template <typename AlphabetType, typename ContextType>
void pps(AlphabetType const* text,
         const uint64_t size,
         const uint64_t levels,
         ContextType& ctx) {
  auto sorted_text_ = std::vector<AlphabetType>(size);
  auto sorted_text = span<AlphabetType>(sorted_text_);

  std::vector<uint64_t> offsets_(1 << levels, 0);
  auto offsets = span<uint64_t>(offsets_);
  auto& bv = ctx.bv();
  auto&& zeros = ctx.zeros();

  #pragma omp parallel
  {
    const auto omp_rank = omp_get_thread_num();
    const auto omp_size = omp_get_num_threads();
    const uint64_t alphabet_size = (1 << levels);

    {
      auto&& rank_hist = ctx.hist_at_shard(omp_rank);

      // While initializing the histogram, we also compute the first level
      omp_write_bits_wordwise(0, size, bv[0], [&](uint64_t const i) {
        rank_hist[text[i]]++;
        uint64_t bit = ((text[i] >> (levels - 1)) & 1ULL);
        return bit;
      });
    }

    // The number of 0's at the last level is the number of "even" characters
    #pragma omp single
    for (int32_t rank = 0; rank < omp_size; ++rank) {
      auto&& rank_hist = ctx.hist_at_shard(rank);
      for (uint64_t i = 0; i < alphabet_size; i += 2) {
        zeros[levels - 1] += rank_hist[i];
      }
    }

    // Now we compute the wavelet structure bottom-up, i.e., the last level
    // first
    for (uint64_t level = levels - 1; level > 0; --level) {
      const uint64_t prefix_shift = (levels - level);
      const uint64_t cur_bit_shift = prefix_shift - 1;

      // Compute the histogram and the border for each bit prefix and
      // processor, i.e., for one fixed bit prefix we compute the prefix sum
      // over the number of occurrences at each processor
      #pragma omp for
      for (uint64_t i = 0; i < alphabet_size; i += (1ULL << prefix_shift)) {
        {
          auto&& hist = ctx.hist_at_shard(0);
          auto&& borders = ctx.borders_at_shard(0);

          hist[i] += hist[i + (1ULL << cur_bit_shift)];
          borders[i] = 0;
        }
        for (int32_t rank = 1; rank < omp_size; ++rank) {
          auto&& hist = ctx.hist_at_shard(rank);
          auto&& borders = ctx.borders_at_shard(rank);
          auto&& prev_borders = ctx.borders_at_shard(rank - 1);
          auto&& prev_hist = ctx.hist_at_shard(rank - 1);

          hist[i] += hist[i + (1ULL << cur_bit_shift)];
          borders[i] = prev_borders[i] + prev_hist[i];
        }
      }

      // Now we compute the offset for each bit prefix, i.e., the number of
      // lexicographically smaller characters
      #pragma omp single
      {
        auto&& last_hist = ctx.hist_at_shard(omp_size - 1);
        auto&& last_borders = ctx.borders_at_shard(omp_size - 1);
        for (uint64_t i = 1; i < (1ULL << level); ++i) {
          const auto rho = ctx.rho(level, i);
          const auto prev_rho = ctx.rho(level, i - 1);

          offsets[rho << prefix_shift] =
              offsets[prev_rho << prefix_shift] +
              last_borders[prev_rho << prefix_shift] +
              last_hist[prev_rho << prefix_shift];
          if (ContextType::compute_rho) {
            ctx.set_rho(level, i - 1, prev_rho >> 1);
          }
        }
        // The number of 0s is the position of the first 1 at the first
        // processor
        if constexpr (ContextType::compute_zeros) {
          zeros[level - 1] = offsets[1ULL << prefix_shift];
        }
      }
      // We add the offset to the borders (for performance)
      #pragma omp for
      for (int32_t rank = 0; rank < omp_size; ++rank) {
        auto&& borders = ctx.borders_at_shard(rank);
        for (uint64_t i = 0; i < alphabet_size; i += (1ULL << prefix_shift)) {
          borders[i] += offsets[i];
        }
      }

      // We align the borders (in memory) to increase performance by reducing
      // the number of cache misses
      std::vector<uint64_t> borders_aligned_(1ULL << level, 0);
      span<uint64_t> borders_aligned(borders_aligned_);
      {
        auto&& borders = ctx.borders_at_shard(omp_rank);
        for (uint64_t i = 0; i < alphabet_size; i += (1ULL << prefix_shift)) {
          borders_aligned[i >> prefix_shift] = borders[i];
        }
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

      omp_write_bits_wordwise(0, size, bv[level], [&](uint64_t pos) {
        return sorted_text[pos] & 1ULL;
      });
    }
  }
}

/******************************************************************************/
