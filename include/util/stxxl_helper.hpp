/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <stxxl/vector>
#include "util/type_for_bytes.hpp"

template <typename value_type>
using stxxlvector = typename stxxl::VECTOR_GENERATOR<value_type>::result;

template <typename value_type>
using stxxlreader = typename stxxl::VECTOR_GENERATOR<value_type>::result::bufreader_type;

template <typename value_type>
using stxxlwriter = typename stxxl::VECTOR_GENERATOR<value_type>::result::bufwriter_type;

template <int word_width>
using external_vector = stxxlvector<typename type_for_bytes<word_width>::type>;

// simple accessor to simulate 2D array on 1D stxxlvector
// [for testing only, as [] operator is very expensive]
template <typename value_type>
class stxxlvector_offset {

private:
  const stxxlvector<value_type> &vec;
  const uint64_t off;

public:
  stxxlvector_offset(stxxlvector<value_type> &vector, uint64_t offset) 
  : vec(vector), off(offset) {}
  
  inline const value_type operator [](const uint64_t index) const {
    return vec[index + off];
  }
};

/******************************************************************************/
