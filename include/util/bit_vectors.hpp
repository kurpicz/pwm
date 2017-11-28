/*******************************************************************************
 * include/util/bit_vectors.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/common.hpp"
#include "util/flat_two_dim_array.hpp"
#include "util/external_flat_two_dim_array.hpp"

struct bit_vector_sizes {
  static uint64_t level_size(const uint64_t, const uint64_t size) {
    return word_size(size);
  };
}; // struct bit_vector_sizes

using bit_vectors = flat_two_dim_array<uint64_t, bit_vector_sizes>;
using external_bit_vectors =
  external_flat_two_dim_array<uint64_t, bit_vector_sizes>;

template <bool is_tree>
struct bit_vectors_switch { };

template <>
struct bit_vectors_switch<true> {
  using type = bit_vectors;
};

template <>
struct bit_vectors_switch<false> {
  using type = external_bit_vectors;
};

template <bool is_tree>
using bit_vectors_t = typename bit_vectors_switch<is_tree>::type;

/******************************************************************************/
