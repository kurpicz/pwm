/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include <type_traits>

#include "util/bit_vectors.hpp"

class wavelet_structure {

public:
  inline wavelet_structure() = default;

  inline wavelet_structure(internal_bit_vectors&& bvs, std::vector<uint64_t>&& zeros):
    bvs_(std::move(bvs)), zeros_(std::move(zeros)) {}

  inline wavelet_structure(internal_bit_vectors&& bvs) :
    bvs_(std::move(bvs)) {}

  // Prevent accidental copies
  inline wavelet_structure(wavelet_structure const& other) = delete;
  inline wavelet_structure& operator =(wavelet_structure const& other) = delete;

  // Allow moving
  inline wavelet_structure(wavelet_structure&& other) = default;
  inline wavelet_structure& operator =(wavelet_structure&& other) = default;

  inline uint64_t levels() const {
    return bvs_.levels();
  }

  inline internal_bit_vectors const& bvs() const {
    return bvs_;
  }

  inline std::vector<uint64_t> const& zeros() const {
    return zeros_;
  }
private:
  internal_bit_vectors bvs_;
  std::vector<uint64_t> zeros_;
}; // class wavelet_structure

/******************************************************************************/
