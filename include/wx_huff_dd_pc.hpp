/*******************************************************************************
 * include/wx_huff_dd_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "huffman/huff_dd.hpp"
#include "huffman/huff_pc.hpp"

#include "wx_huff_pc.hpp"

struct huff_pc_disp {
  static constexpr bool needs_all_borders = false;

  template <typename AlphabetType, typename ContextType, typename HuffCodes>
  static void calc_huff(AlphabetType const* text,
                        uint64_t const size,
                        uint64_t const levels,
                        HuffCodes const& codes,
                        ContextType& ctx,
                        span<uint64_t const> const level_sizes) {
    huff_pc(text, size, levels, codes, ctx, level_sizes);
  }
};

template <typename AlphabetType, bool is_tree_>
using wx_huff_dd_pc = huff_dd<huff_pc_disp, wx_huff_pc, AlphabetType, is_tree_>;

/******************************************************************************/
