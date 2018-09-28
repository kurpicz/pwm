/*******************************************************************************
 * include/huffman/huff_dd.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <cstring>
#include <omp.h>
#include <vector>

#include "construction/wavelet_structure.hpp"
#include "arrays/slice.hpp"

#include "huffman/huff_merge.hpp"
#include "huffman/huff_bit_vectors.hpp"
#include "huffman/huff_codes.hpp"
#include "huffman/ctx_huff_all_levels.hpp"
#include "huffman/huff_parallel_level_sizes_builder.hpp"

template <typename Algorithm, typename AlphabetType, bool is_tree_>
class huff_dd {

public:
  static constexpr bool  is_parallel = true;
  static constexpr bool  is_tree   = is_tree_;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);
  static constexpr bool  is_huffman_shaped = true;

  using ctx_t = ctx_huff_all_levels<is_tree>;

  template <typename InputType>
  static wavelet_structure compute(const InputType& global_text,
    const uint64_t size, const uint64_t /*levels*/) {

    // TODO ^ remove levels parameter from API

    const uint64_t shards = omp_get_max_threads();

    std::vector<histogram<AlphabetType>> local_hists(shards);

    auto get_local_text = [global_text, size, shards] (size_t rank) {
      const uint64_t local_size = (size / shards) +
        ((rank < size % shards) ? 1 : 0);
      const uint64_t offset = (rank * (size / shards)) +
        std::min<uint64_t>(rank, size % shards);

      const AlphabetType* local_text = global_text + offset;

      return Slice<AlphabetType const> { local_text, local_size };
    };

    #pragma omp parallel for
    for (size_t shard = 0; shard < shards; shard++) {
      auto text = get_local_text(shard);

      // calculate local histogram
      local_hists[shard] = histogram<AlphabetType>(text.data(), text.size());
    }

    // merge local histograms
    histogram<AlphabetType> merged_hist(local_hists);

    // prepare level size builder
    auto builder = parallel_level_sizes_builder<AlphabetType> {
        std::move(merged_hist),
        std::move(local_hists),
    };

    // build huffman codes and calculate level sizes
    canonical_huff_codes<AlphabetType, is_tree> codes
      = canonical_huff_codes<AlphabetType, is_tree>(builder);
    uint64_t const levels = builder.levels();

    if(size == 0) {
      if constexpr (ctx_t::compute_zeros) {
        return wavelet_structure_matrix_huffman(std::move(codes));
      } else {
        return wavelet_structure_tree_huffman(std::move(codes));
      }
    }

    const auto rho = rho_dispatch<is_tree>::create(levels);
    auto ctxs = std::vector<ctx_t>(shards);

    #pragma omp parallel for
    for (size_t shard = 0; shard < shards; shard++) {
      auto const text = get_local_text(shard);
      std::vector<std::vector<uint64_t>> const& local_level_sizes
        = builder.level_sizes_shards();

      ctxs[shard] = ctx_t(local_level_sizes[shard], levels, rho);

      Algorithm::calc_huff(text.data(), text.size(), levels, codes, ctxs[shard]);

      ctxs[shard].discard_non_merge_data();
    }

    std::vector<uint64_t> const& level_sizes = builder.level_sizes();
    auto bv = huff_merge_bit_vectors(level_sizes, shards, ctxs, rho);

    if constexpr (ctx_t::compute_zeros) {
      auto zeros = std::vector<uint64_t>(levels, 0);

      for(size_t level = 0; level < levels; level++) {
        for(size_t shard = 0; shard < ctxs.size(); shard++) {
          zeros[level] += ctxs[shard].zeros()[level];
        }
      }

      return wavelet_structure_matrix_huffman(
        std::move(bv), std::move(zeros), std::move(codes));
    } else {
      return wavelet_structure_tree_huffman(
        std::move(bv), std::move(codes));
    }
  }
};

/******************************************************************************/
