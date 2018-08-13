/*******************************************************************************
 * include/util/memory_types.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/common.hpp"
#include "util/bit_vectors.hpp"
#include "util/flat_two_dim_array.hpp"
#include "util/flat_two_dim_array_external.hpp"
#include "util/type_for_bytes.hpp"
#include "util/stxxl_helper.hpp"
#include "util/wavelet_structure.hpp"


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
  using type = wavelet_structure;
};

template <>
struct out_type<true> {
  using type = external_bit_vectors;
};

enum memory_mode { internal, external_input, external_output, external };

//~ template <memory_mode mem_mode, int bytes_per_word>
//~ struct type_info {
  //~ static constexpr bool input_external = 
    //~ mem_mode == memory_mode::external || 
    //~ mem_mode == memory_mode::external_input;
    
  //~ static constexpr bool output_external = 
    //~ mem_mode == memory_mode::external || 
    //~ mem_mode == memory_mode::external_output;

  //~ using input_type = in_type<input_external, bytes_per_word>;
  //~ using output_type = out_type<output_external>;
//~ };

/******************************************************************************/
