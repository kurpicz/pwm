/*******************************************************************************
 * include/wx_huff_naive.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>

#include "construction/wavelet_structure.hpp"

#include "huffman/ctx_huff_all_levels.hpp"
#include "huffman/huff_bit_vectors.hpp"
#include "huffman/huff_codes.hpp"
#include "huffman/huff_level_sizes_builder.hpp"
#include "huffman/huff_naive.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_huff_naive {

public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = true;

  // TODO: change to single level
  using ctx_t = ctx_huff_all_levels<is_tree>;

  static wavelet_structure compute(AlphabetType const* const text,
                                   const uint64_t size,
                                   const uint64_t /*levels*/) {

    histogram<AlphabetType> hist{text, size};
    level_sizes_builder<AlphabetType> builder{std::move(hist)};
    canonical_huff_codes<AlphabetType, is_tree> codes(builder);

    auto const& level_sizes = builder.level_sizes();
    uint64_t const levels = builder.levels();

    if (size == 0) {
      if constexpr (is_tree) {
        return wavelet_structure_tree_huffman<AlphabetType>(std::move(codes));
      } else {
        return wavelet_structure_matrix_huffman<AlphabetType>(std::move(codes));
      }
    }

    const auto rho = rho_dispatch<is_tree>::create(levels);
    auto ctx = ctx_t(level_sizes, levels, rho);

    huff_naive(text, size, levels, codes, ctx);

    auto& bv = ctx.bv();
    auto& zeros = ctx.zeros();

    if constexpr (is_tree) {
      return wavelet_structure_tree_huffman<AlphabetType>(std::move(bv),
                                                          std::move(codes));
    } else /*if constexpr (!is_tree)*/ {
      return wavelet_structure_matrix_huffman<AlphabetType>(
          std::move(bv), std::move(zeros), std::move(codes));
    }
  }
}; // class wx_huff_naive

/******************************************************************************/
