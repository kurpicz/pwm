/*******************************************************************************
 * include/huffman/huff_border_hist_array.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/common.hpp"
#include "arrays/flat_two_dim_array.hpp"

struct huff_border_hist_level_size {
  static uint64_t level_size(const uint64_t level) {
    return 1ULL << level;
  };

  static constexpr bool is_bit_vector = false;
}; // struct huff_border_hist_level_size

using huff_border_hist_array = flat_two_dim_array<uint64_t, huff_border_hist_level_size>;

/******************************************************************************/
