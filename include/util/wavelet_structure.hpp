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
#include "huffman/huff_bit_vectors.hpp"

class wavelet_structure {

public:
  inline wavelet_structure() = default;

  inline wavelet_structure(bit_vectors&& bvs, std::vector<uint64_t>&& zeros)
  : bvs_(std::move(bvs)), zeros_(std::move(zeros)) { }

  wavelet_structure(huff_bit_vectors&& bvs, std::vector<uint64_t>&& zeros)
  : huff_bvs_(std::move(bvs)), zeros_(std::move(zeros)) { }

  inline wavelet_structure(bit_vectors&& bvs) : bvs_(std::move(bvs)) {}

  wavelet_structure(huff_bit_vectors&& bvs) : huff_bvs_(std::move(bvs)) { }

  // Prevent accidental copies
  wavelet_structure(wavelet_structure const&) = delete;
  wavelet_structure& operator =(wavelet_structure const&) = delete;

  // Allow moving
  inline wavelet_structure(wavelet_structure&& other) = default;
  inline wavelet_structure& operator =(wavelet_structure&& other) = default;

  inline uint64_t levels() const {
    return bvs_.levels();
  }

  inline bit_vectors const& bvs() const {
    return bvs_;
  }

  huff_bit_vectors const& huff_bvs() const {
    return huff_bvs_;
  }

  inline std::vector<uint64_t> const& zeros() const {
    return zeros_;
  }
private:
  bit_vectors bvs_;
  huff_bit_vectors huff_bvs_;
  std::vector<uint64_t> zeros_;
}; // class wavelet_structure

/******************************************************************************/
