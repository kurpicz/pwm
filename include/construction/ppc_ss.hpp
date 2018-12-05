/*******************************************************************************
 * include/construction/ppc_ss.hpp
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
void ppc_ss(AlphabetType const* text,
            const uint64_t size,
            const uint64_t levels,
            ContextType& ctx) {
  auto& bv = ctx.bv();
  auto&& zeros = ctx.zeros();

  const uint64_t alphabet_size = (1 << levels);

  std::vector<std::vector<uint64_t>> all_hists_ =
      std::vector<std::vector<uint64_t>>(omp_get_max_threads());
  span<std::vector<uint64_t>> all_hists(all_hists_);

  #pragma omp parallel
  {
    const auto omp_rank = omp_get_thread_num();
    all_hists[omp_rank] = std::vector<uint64_t>(alphabet_size + 1, 0);
    auto&& hist = span<uint64_t>(all_hists[omp_rank]);

    omp_write_bits_wordwise(0, size, bv[0], [&](uint64_t i) {
      hist[text[i]]++;
      uint64_t const bit = ((text[i] >> (levels - 1)) & 1ULL);
      return bit;
    });
  }

  auto&& hist = ctx.hist_at_level(levels);
  for (uint64_t i = 0; i < all_hists.size(); ++i) {
    for (uint64_t j = 0; j < alphabet_size; ++j) {
      hist[j] += all_hists[i][j];
    }
  }
  if constexpr (ContextType::compute_zeros) {
    for (uint64_t i = 0; i < alphabet_size; i += 2) {
      zeros[levels - 1] += hist[i];
    }
  }

  bottom_up_compute_hist_and_borders_and_optional_zeros(size, levels, ctx);

  #pragma omp parallel num_threads(levels)
  {
    uint64_t level = omp_get_thread_num();
    auto&& borders = ctx.borders_at_level(level);
    for (uint64_t i = 0; i < size; ++i) {
      single_scan_write_bit(bv, level, levels, borders, text[i]);
    }
  }
}

/******************************************************************************/
