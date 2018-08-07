/*******************************************************************************
 * include/wx_dd_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <cstring>
#include <omp.h>
#include <vector>

#include "util/ctx_all_levels.hpp"
#include "util/merge.hpp"
#include "util/pc.hpp"
#include "util/wavelet_structure.hpp"
#include "util/memory_types.hpp"

template <typename AlphabetType, bool is_tree_, memory_mode mem_mode_>
class wx_dd_pc {

public:
  static constexpr bool  is_parallel = true;
  static constexpr bool  is_tree   = is_tree_;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);
  static constexpr bool  is_huffman_shaped = false;
  static constexpr memory_mode mem_mode = mem_mode_;

  template <typename InputType, typename OutputType>
  static wavelet_structure<OutputType> compute(const InputType& global_text,
    const uint64_t size, const uint64_t levels) {
        
    using ctx_t = ctx_all_levels<OutputType, is_tree>;
    
    if(size == 0) { return wavelet_structure<OutputType>(); }

    const uint64_t shards = omp_get_max_threads();

    const auto rho = rho_dispatch<is_tree>::create(levels);
    auto ctxs = std::vector<ctx_t>(shards);

    for (size_t shard = 0; shard < shards; shard++) {
      const uint64_t local_size = (size / shards) +
        ((shard < size % shards) ? 1 : 0);

      // TODO: store size in ctx so that later can reload it
      ctxs[shard] = ctx_t(local_size, levels, rho);
    }

    #pragma omp parallel
    {
      const uint64_t omp_rank = omp_get_thread_num();
      const uint64_t omp_size = omp_get_num_threads();
      assert(omp_size == shards);

      const uint64_t local_size = (size / omp_size) +
        ((omp_rank < size % omp_size) ? 1 : 0);
      const uint64_t offset = (omp_rank * (size / omp_size)) +
        std::min<uint64_t>(omp_rank, size % omp_size);

      const InputType text = global_text + offset;

      pc(text, local_size, levels, ctxs[omp_rank]);
    }

    for(auto& ctx : ctxs) {
      ctx.discard_non_merge_data();
    }
    auto _bv = merge_bit_vectors(size, levels, shards, ctxs, rho);

    if (ctx_t::compute_zeros) {
      auto _zeros = std::vector<uint64_t>(levels, 0);

      #pragma omp parallel for
      for(size_t level = 0; level < levels; level++) {
        for(size_t shard = 0; shard < ctxs.size(); shard++) {
          _zeros[level] += ctxs[shard].zeros()[level];
        }
      }

      return wavelet_structure<OutputType>(std::move(_bv), std::move(_zeros));
    } else {
      return wavelet_structure<OutputType>(std::move(_bv));
    }
  }
}; // class wx_dd_pc

/******************************************************************************/
