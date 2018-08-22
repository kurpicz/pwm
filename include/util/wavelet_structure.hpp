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

#include "flat_two_dim_array.hpp"

using base_bit_vectors = base_flat_two_dim_array<uint64_t>;

class wavelet_structure {

public:
  wavelet_structure() = default;

  wavelet_structure(base_bit_vectors&& bvs, std::vector<uint64_t>&& zeros)
  : bvs_(std::move(bvs)), zeros_(std::move(zeros)) { }

  wavelet_structure(base_bit_vectors&& bvs) : bvs_(std::move(bvs)) { }

  // Prevent accidental copies
  wavelet_structure(wavelet_structure const&) = delete;
  wavelet_structure& operator =(wavelet_structure const&) = delete;

  // Allow moving
  wavelet_structure(wavelet_structure&& other) = default;
  wavelet_structure& operator =(wavelet_structure&& other) = default;

  inline uint64_t levels() const {
    return bvs_.levels();
  }

  inline base_bit_vectors const& bvs() const {
    return bvs_;
  }

  inline std::vector<uint64_t> const& zeros() const {
    return zeros_;
  }
private:
  base_bit_vectors bvs_;
  std::vector<uint64_t> zeros_;
}; // class wavelet_structure

/******************************************************************************/
