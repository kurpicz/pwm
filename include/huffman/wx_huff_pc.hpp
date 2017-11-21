/*******************************************************************************
 * include/huffman/wx_huff_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>

#include "ctx_huff_all_levels.hpp"
#include "huff_codes.hpp"
#include "util/alphabet_util.hpp"
#include "util/histogram.hpp"

template <typename AlphabetType, bool is_matrix>
class wx_huff_pc {
  using ctx_t = ctx_huff_all_levels<is_matrix>;

public:
  static constexpr bool is_parallel = false;
  static constexpr bool is_tree = !is_matrix;
  static constexpr uint8_t word_width = sizeof(AlphabetType);

  static wavelet_structure compute(AlphabetType const* const text,
    const uint64_t size, const uint64_t reduced_sigma) {

    if (size == 0) {
      return wavelet_structure();
    }

    std::vector<uint64_t> hist =
      compute_initial_histogram(text, size, reduced_sigma);
    auto codes = huff_codes<AlphabetType, is_matrix>(hist);
    // There is no code with length longer than the number of codes 
    std::vector<uint64_t> sizes(codes.code_pairs().size(), 0);

    // TODO: This can be done during code computation
    uint64_t max_code_length = 0;
    for (uint64_t i = 0; i < size; ++i) {
      const auto code_length = codes.encode_symbol(text[i]).code_length;
      max_code_length = std::max(max_code_length, code.code_length);
      ++sizes[code_length];
    }
    sizes.resize(max_code_length);
    for (uint64_t i = sizes.size() - 1; i > 0; --i) {
      sizes[i - 1] += sizes[i];
    }
    assert(sizes[0] == size);

    const uint64_t levels = levels_for_max_char(sizes.back());
    const auto rho = rho_dispatch<is_matrix>::create(levels);

    auto ctx = ctx_t(sizes, levels, rho)
    huff_pc(text, size, ctx, codes);

    if (ctx_t::compute_zeros) {
      return wavelet_structure(std::move(ctx.bv()), std::move(ctx.zeros()));
    } else {
      return wavelet_structure(std::move(ctx.bv()));
    }
  }
}; // class wx_huff_pc

/******************************************************************************/
