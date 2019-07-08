
/*******************************************************************************
 * include/wx_pps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "construction/building_blocks.hpp"
#include "construction/ctx_generic.hpp"
#include "construction/pps.hpp"
#include "construction/wavelet_structure.hpp"
#include "util/common.hpp"

#include "wx_base.hpp"
#include "wx_ps.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_pps : public wx_in_out_external<false, false>  {

public:
  static constexpr bool is_parallel = true;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  using ctx_t = ctx_generic<is_tree,
                            ctx_options::borders::sharded_single_level,
                            ctx_options::hist::sharded_single_level,
                            ctx_options::live_computed_rho,
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

    const uint64_t shards = omp_get_max_threads();

    if (shards <= 1) {
      return wx_ps<AlphabetType, is_tree>::
        template compute<InputType>(text,
                                    size,
                                    levels);
    }

    ctx_t ctx(size, levels, levels, shards);

    pps(text, size, levels, ctx);

    if constexpr (ctx_t::compute_zeros) {
      return wavelet_structure_matrix(std::move(ctx.bv()),
                                      std::move(ctx.take_zeros()));
    } else {
      return wavelet_structure_tree(std::move(ctx.bv()));
    }
  }
};

/******************************************************************************/
