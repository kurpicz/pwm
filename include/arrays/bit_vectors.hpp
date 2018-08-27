/*******************************************************************************
 * include/arrays/bit_vectors.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "flat_two_dim_array.hpp"

struct bit_vector_sizes {
  static uint64_t level_size(const uint64_t, const uint64_t size) {
    return size;
  };

  static constexpr bool is_bit_vector = true;
}; // struct bit_vector_sizes

struct huff_bit_vectors_sizes {
  static uint64_t level_size(const uint64_t level,
    std::vector<uint64_t> const& level_sizes) {
    return level_sizes[level];
  }

  static constexpr bool is_bit_vector = true;
}; // struct huff_bit_vectors_sizes

using base_bit_vectors = base_flat_two_dim_array<uint64_t>;
using bit_vectors = flat_two_dim_array<uint64_t, bit_vector_sizes>;
using huff_bit_vectors = flat_two_dim_array<uint64_t, huff_bit_vectors_sizes>;

/******************************************************************************/
