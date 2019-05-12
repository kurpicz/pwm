/*******************************************************************************
 * include/util/memory_types.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/type_for_bytes.hpp"

#include "external_memory/stxxl_helper.hpp"

#include "construction/wavelet_structure.hpp"
#include "external_memory/wavelet_structure_external.hpp"

enum memory_mode : bool { internal = false, external = true };

template <bool external, int word_width>
struct input_type;

template <bool external>
struct output_type;

template <int word_width>
struct input_type<memory_mode::internal, word_width> {
  using value_type = typename type_for_bytes<word_width>::type;
  using type = value_type*;
};

template <int word_width>
struct input_type<memory_mode::external, word_width> {
  using value_type = typename type_for_bytes<word_width>::type;
  using type = stxxlvector<value_type>;
};

template <>
struct output_type<memory_mode::internal> {
  using type = wavelet_structure;
};

template <>
struct output_type<memory_mode::external> {
  using type = wavelet_structure_external;
};

/******************************************************************************/
