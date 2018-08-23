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

  wavelet_structure(base_bit_vectors&& bvs,
                    std::vector<uint64_t>&& zeros,
                    bool is_huffman_shaped)
  : bvs_(std::move(bvs))
  , zeros_(std::move(zeros))
  , is_tree_(false)
  , is_huffman_shaped_(is_huffman_shaped) { }

  wavelet_structure(base_bit_vectors&& bvs,
                    bool is_huffman_shaped)
  : bvs_(std::move(bvs))
  , is_tree_(true)
  , is_huffman_shaped_(is_huffman_shaped) { }

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

  inline bool is_tree() const { return is_tree_; }
  inline bool is_huffman_shaped() const { return is_huffman_shaped_; }
private:
  base_bit_vectors bvs_;
  std::vector<uint64_t> zeros_;
  bool is_tree_;
  bool is_huffman_shaped_;
}; // class wavelet_structure

/******************************************************************************/
