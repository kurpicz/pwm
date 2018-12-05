/*******************************************************************************
 * include/construction/ppc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <omp.h>

#include "arrays/span.hpp"
#include "construction/building_blocks.hpp"
#include "construction/wavelet_structure.hpp"
#include "util/common.hpp"

template <typename AlphabetType, typename ContextType>
void ppc(AlphabetType const* text,
         const uint64_t size,
         const uint64_t levels,
         ContextType& ctx) {
  auto&& bv = ctx.bv();

  // compute initial histogram
  {
    const uint64_t alphabet_size = (1ull << levels);
    std::vector<uint64_t> initial_hist_v(alphabet_size * levels, 0);
    span<uint64_t> initial_hist{initial_hist_v.data(), initial_hist_v.size()};
    auto&& hist = ctx.hist_at_level(levels);

    #pragma omp parallel num_threads(levels)
    {
      const uint64_t omp_rank = uint64_t(omp_get_thread_num());
      DCHECK(levels == uint64_t(omp_get_num_threads()));

      auto initial_hist_ptr = initial_hist.slice(
          alphabet_size * omp_rank, alphabet_size * (omp_rank + 1));

      omp_write_bits_wordwise(0, size, bv[0], [&](uint64_t const i) {
        auto const c = text[i];
        initial_hist_ptr[c]++;
        uint64_t bit = ((c >> (levels - 1)) & 1ULL);
        return bit;
      });

      #pragma omp barrier

      // Compute the histogram with respect to the local slices of the text
      #pragma omp for
      for (uint64_t i = 0; i < alphabet_size; ++i) {
        for (uint64_t rank = 0; rank < levels; ++rank) {
          hist[i] += initial_hist[(alphabet_size * rank) + i];
        }
      }
    }
  }

  auto&& leaf_hist = ctx.hist_at_level(levels);

  #pragma omp parallel num_threads(levels)
  {
    const uint64_t omp_rank_ = omp_get_thread_num();
    DCHECK(omp_rank_ < levels);

    // Compute the histogram for each level of the wavelet structure
    #pragma omp for
    for (uint64_t level = 1; level < levels; ++level) {
      // TODO: remove again
      const uint64_t omp_rank = omp_get_thread_num();
      DCHECK(omp_rank < levels);
      DCHECK(omp_rank == omp_rank_);

      auto&& hist = ctx.hist_at_level(level);
      const uint64_t blocks = (1 << level);
      const uint64_t required_characters = (1 << (levels - level));
      for (uint64_t i = 0; i < blocks; ++i) {
        for (uint64_t j = 0; j < required_characters; ++j) {
          hist[i] += leaf_hist[(i * required_characters) + j];
        }
      }
    }

    // Now we compute the wavelet structure bottom-up, i.e., the last level
    // first
    #pragma omp for
    for (uint64_t level = 1; level < levels; ++level) {
      // TODO: remove again
      const uint64_t omp_rank = omp_get_thread_num();
      DCHECK(omp_rank < levels);
      DCHECK(omp_rank == omp_rank_);

      const uint64_t blocks = (1 << level);
      const uint64_t prefix_shift = (levels - level);
      const uint64_t cur_bit_shift = prefix_shift - 1;

      #pragma omp critical
      {
        auto borders = ctx.borders_at_shard(omp_rank).slice(0, blocks);

        compute_borders_and_optional_zeros_and_optional_rho(level, blocks, ctx,
                                                            borders);

        // Now we insert the bits with respect to their bit prefixes

        for (uint64_t i = 0; i < size; ++i) {
          const uint64_t pos = borders[text[i] >> prefix_shift]++;
          bv[level][pos >> 6] |=
              (((text[i] >> cur_bit_shift) & 1ULL) << (63ULL - (pos & 63ULL)));
        }
      }
    }
  }

  ctx.hist_at_level(0)[0] = size;
  ctx.discard_borders();
}

/******************************************************************************/
