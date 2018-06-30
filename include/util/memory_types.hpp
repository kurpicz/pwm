/*******************************************************************************
 * include/util/bit_vectors.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <stxxl/vector>

#include "util/common.hpp"
#include "util/flat_two_dim_array.hpp"
#include "util/external_flat_two_dim_array.hpp"
#include "util/wavelet_structure.hpp"
#include "util/type_for_bytes.hpp"
#include "util/stxxl_helper.hpp"

using internal_bit_vectors = flat_two_dim_array<uint64_t, bit_vector_sizes>;
using external_bit_vectors =
  external_flat_two_dim_array<uint64_t, bit_vector_sizes>;

template <bool external, int bytes_per_word>
struct in_type { };

template <bool external>
struct out_type { };

template <int bytes_per_word>
struct in_type<false, bytes_per_word> {
  using type = typename type_for_bytes<bytes_per_word>::type const *;
};

template <int bytes_per_word>
struct in_type<true, bytes_per_word> {
  using type = stxxlvector<typename type_for_bytes<bytes_per_word>::type> const;
};

template <>
struct out_type<false> {
  using type = internal_bit_vectors;
};

template <>
struct out_type<true> {
  using type = external_bit_vectors;
};


/******************************************************************************/
