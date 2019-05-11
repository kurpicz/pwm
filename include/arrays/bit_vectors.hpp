/*******************************************************************************
 * include/arrays/bit_vectors.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "flat_two_dim_array.hpp"

template <bool requires_initialization_>
struct bit_vectors_config {
  static uint64_t level_size(const uint64_t, const uint64_t size) {
    return size;
  };

  static constexpr bool is_bit_vector = true;
  static constexpr bool requires_initialization = requires_initialization_;
}; // struct bit_vectors_config

using base_bit_vectors = base_flat_two_dim_array<uint64_t>;

template <bool requires_initialization = true>
using bit_vectors =
    flat_two_dim_array<uint64_t, bit_vectors_config<requires_initialization>>;

/******************************************************************************/
