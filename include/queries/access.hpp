/*******************************************************************************
 * include/queries/access.hpp
 *
 * Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>

#include "construction/wavelet_structure.hpp"
#include "queries/bin_rank_popcnt.hpp"
#include "util/common.hpp"

template <typename rank_support_type = bin_rank_popcnt>
class access_support {

public:
  access_support(wavelet_structure& ws) : ws_(ws), rank_support_(ws_.levels()) {
    for (size_t i = 0; i < ws_.levels(); ++i) {
      rank_support_[i] = std::move(rank_support_type(ws_.bvs()[i],
        word_size(ws_.bvs().level_bit_size(i))));
    }
  }

  // Returns the bit pattern of the character that doesn't necessarily
  // correspond to the encoding of the character, as we use an effective
  // alphabet.
  // TODO: make the effective alphabet accessible for this, such that we can
  // transform the output!
  inline uint64_t operator [](size_t index) const {
    if (ws_.is_tree() && !ws_.is_huffman_shaped()) {
      return access_tree(0, index, 0, ws_.bvs().level_bit_size(0), 0ULL) ;}
    else if (!ws_.is_tree() && !ws_.is_huffman_shaped()) {
      return access_matrix(0, index, 0ULL);
    } else { std::cout << "NOT YET IMPLEMENTED!" << std::endl; }
    return 0;
  }

private:
  inline uint64_t access_tree(const size_t level, const size_t index,
    const size_t start, const size_t end, uint64_t word) const {

    if (level == ws_.levels()) { return word; }

    size_t const left = rank_support_[level].rank0(start);
    size_t const right = rank_support_[level].rank0(end);

    word <<= 1;
    if (bit_at(ws_.bvs()[level], start + index)) {
      word |= 1ULL;
      size_t const offset = rank_support_[level].rank1(start + index);
      return access_tree(level + 1, offset - (start - left),
        start + right - left, end, word);
    } else {
      size_t const offset = rank_support_[level].rank0(start + index);
      return access_tree(level + 1, offset - left, start, start + right - left,
        word);
    }
    
  }

  inline uint64_t access_matrix(size_t const level, size_t index,
    uint64_t word) const {

    if (level == ws_.levels()) { return word; }

    word <<= 1;
    if (bit_at(ws_.bvs()[level], index)) {
      word |= 1ULL;
      index = rank_support_[level].rank1(index) + ws_.zeros()[level];
    } else { index = rank_support_[level].rank0(index); }
    return access_matrix(level + 1, index, word);
  }

private:
  wavelet_structure& ws_;
  std::vector<rank_support_type> rank_support_;  

}; // class access_support



/******************************************************************************/
