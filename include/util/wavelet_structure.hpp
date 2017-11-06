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
#include "bit_vectors.hpp"

class wavelet_structure {

public:
  inline wavelet_structure() = default;

  inline wavelet_structure(bit_vectors&& bvs, std::vector<uint64_t>&& zeros):
    bvs_(std::move(bvs)), zeros_(std::move(zeros)) {}

  inline wavelet_structure(bit_vectors&& bvs) :
    bvs_(std::move(bvs)) {}

  // Prevent accidental copies
  inline wavelet_structure(wavelet_structure const& other) = delete;
  inline wavelet_structure& operator =(wavelet_structure const& other) = delete;

  // Allow moving
  inline wavelet_structure(wavelet_structure&& other) = default;
  inline wavelet_structure& operator =(wavelet_structure&& other) = default;

  inline std::vector<uint64_t*> const& raw_bvs() const {
    return bvs_.vec();
  }

  inline std::vector<uint64_t> const& raw_zeros() const {
    return zeros_;
  }
private:
  bit_vectors bvs_;
  std::vector<uint64_t> zeros_;
}; // class wavelet_structure

/******************************************************************************/
