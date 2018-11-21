/*******************************************************************************
 * include/arrays/helper_array.hpp
 *
 * Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

struct helper_array_config {
  static uint64_t level_size(const uint64_t, const uint64_t size) {
    return size;
  }

  static constexpr bool is_bit_vector = false;
  static constexpr bool requires_initialization = true;
}; // struct helper_array_config

using helper_array = flat_two_dim_array<uint64_t, helper_array_config>;

/******************************************************************************/
