/*******************************************************************************
 * include/wx_huff_dd_naive.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "huffman/huff_dd.hpp"
#include "huffman/huff_naive.hpp"

#include "wx_huff_naive.hpp"

struct huff_naive_disp {
  static constexpr bool needs_all_borders = false;

  template <typename AlphabetType, typename ContextType, typename HuffCodes>
  static void calc_huff(AlphabetType const* text,
                        uint64_t const size,
                        uint64_t const levels,
                        HuffCodes const& codes,
                        ContextType& ctx,
                        span<uint64_t const> const level_sizes) {
    huff_naive(text, size, levels, codes, ctx, level_sizes);
  }
};

template <typename AlphabetType, bool is_tree_>
using wx_huff_dd_naive = huff_dd<huff_naive_disp,
                                 wx_huff_naive,
                                 AlphabetType,
                                 is_tree_>;

/******************************************************************************/
