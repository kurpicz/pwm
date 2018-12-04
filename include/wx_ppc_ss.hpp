/*******************************************************************************
 * include/wx_ppc_ss.hpp
 *
 * Copyright (C) 2017-2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <omp.h>

#include "construction/building_blocks.hpp"
#include "construction/ctx_generic.hpp"
#include "construction/pc_ss.hpp"
#include "construction/wavelet_structure.hpp"

#include "wx_base.hpp"

template <typename AlphabteType, bool is_tree_>
class wx_ppc_ss : public wx_in_out_external<false, false> {

public:
  static constexpr bool is_parallel = true;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabteType);
  static constexpr bool is_huffman_shaped = false;

  using ctx_t = ctx_generic<is_tree,
                            ctx_options::borders::all_level,
                            ctx_options::hist::all_level,
                            ctx_options::pre_computed_rho,
                            ctx_options::bv_initialized,
                            bit_vectors>;

  template <typename InputType>
  static wavelet_structure
  compute(const InputType& text, const uint64_t size, const uint64_t levels) {

    if (size == 0) {
      if constexpr (ctx_t::compute_zeros) {
        return wavelet_structure_matrix();
      } else {
        return wavelet_structure_tree();
      }
    }

    const auto rho = rho_dispatch<is_tree>::create(levels);
    auto ctx = ctx_t(size, levels, rho);
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
    if constexpr (ctx_t::compute_zeros) {
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

    if constexpr (ctx_t::compute_zeros) {
      return wavelet_structure_matrix(std::move(ctx.bv()),
                                      std::move(ctx.take_zeros()));
    } else {
      return wavelet_structure_tree(std::move(ctx.bv()));
    }
  }
}; // class wx_ppc_ss

/******************************************************************************/
