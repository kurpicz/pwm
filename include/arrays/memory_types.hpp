/*******************************************************************************
 * include/util/memory_types.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/type_for_bytes.hpp"

#include "memory_modes.hpp"
#include "stxxl_helper.hpp"

#include "construction/wavelet_structure.hpp"
#include "construction/wavelet_structure_external.hpp"

template <bool external, int word_width>
struct input_type;

template <bool external>
struct output_type;

template <int word_width>
struct input_type<memory_mode::internal, word_width> {
  using type = typename type_for_bytes<word_width>::type*;
};

template <int word_width>
struct input_type<memory_mode::external, word_width> {
  using type = stxxlvector<typename type_for_bytes<word_width>::type>;
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
