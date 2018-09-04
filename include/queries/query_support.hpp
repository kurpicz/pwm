/*******************************************************************************
 * include/queries/query_support.hpp
 *
 * Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>

#include "construction/wavelet_structure.hpp"
#include "queries/bin_rank_popcnt.hpp"
#include "queries/bin_select_popcnt.hpp"
#include "util/common.hpp"

template <typename rank_support_type = bin_rank_popcnt,
          typename select0_support_type = bin_select0_popcnt,
          typename select1_support_type = bin_select1_popcnt>
class query_support {

public:
  query_support(wavelet_structure& ws) : ws_(ws), rank_support_(ws_.levels()) {
    for (size_t i = 0; i < ws_.levels(); ++i) {
      rank_support_[i] = std::move(rank_support_type(ws_.bvs()[i],
        word_size(ws_.bvs().level_bit_size(i))));
      select0_support_.emplace_back(select0_support_type(rank_support_[i]));
      select1_support_.emplace_back(select1_support_type(rank_support_[i]));
    }
  }

  // Returns the bit pattern of the character that doesn't necessarily
  // correspond to the encoding of the character, as we use an effective
  // alphabet.
  // TODO: make the effective alphabet accessible for this, such that we can
  // transform the output!
  inline uint64_t access(size_t index) const {
    if (ws_.is_tree() && !ws_.is_huffman_shaped()) {
      return access_tree(0, index, 0, ws_.bvs().level_bit_size(0), 0ULL);
    }
    else if (!ws_.is_tree() && !ws_.is_huffman_shaped()) {
      return access_matrix(0, index, 0ULL);
    } else { std::cout << "NOT YET IMPLEMENTED!" << std::endl; }
    return 0;
  }

  inline uint64_t rank(uint8_t const symbol, uint64_t const index) const {
    if (index == 0) { return 0; }

    if (ws_.is_tree() && !ws_.is_huffman_shaped()) {
      return rank_tree(0, symbol, index, 0, ws_.bvs().level_bit_size(0));
    }
    else if (!ws_.is_tree() && !ws_.is_huffman_shaped()) {
      return rank_matrix(0, symbol, index, 0);
    } else { std::cout << "NOT YET IMPLEMENTED!" << std::endl; }
    return 0;
  }

  inline uint64_t select(uint8_t const symbol, uint64_t const occ) const {
    if (occ == 0) { return 0; }

    if (ws_.is_tree() && !ws_.is_huffman_shaped()) {
      return select_tree(0, symbol, occ, 0, ws_.bvs().level_bit_size(0));
    }
    else if (!ws_.is_tree() && !ws_.is_huffman_shaped()) {
      return select_matrix(0, symbol, occ, 0);
    } else { std::cout << "NOT YET IMPLEMENTED!" << std::endl; }
    return 0;
  }

private:
  inline uint64_t access_tree(size_t const level, size_t const index,
    size_t const start, size_t const end, uint64_t word) const {

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

  inline uint64_t rank_tree(size_t const level, uint8_t const symbol,
    size_t const index, size_t const start, size_t const end) const {

    if (level == ws_.levels()) { return index; }

    size_t const left = rank_support_[level].rank0(start);
    size_t const right = rank_support_[level].rank0(end);

    size_t const shift_for_bit = ws_.levels() - level - 1;
    if ((symbol >> shift_for_bit) & uint8_t(1)) {
      size_t const offset = rank_support_[level].rank1(start + index);
      return rank_tree(level + 1, symbol, offset - (start - left),
        start + right - left, end);
    } else {
      size_t const offset = rank_support_[level].rank0(start + index);
      return rank_tree(level + 1, symbol, offset - left, start,
        start + right - left);
    }
  }

  inline uint64_t rank_matrix(size_t const level, uint8_t const symbol,
    size_t index, size_t offset) const {

    if (level == ws_.levels()) { return index - offset; }

    size_t const shift_for_bit = ws_.levels() - level - 1;
    if ((symbol >> shift_for_bit) & uint8_t(1)) {
      index = ws_.zeros()[level] + rank_support_[level].rank1(index);
      offset = ws_.zeros()[level] + rank_support_[level].rank1(offset);
    } else {
      index = rank_support_[level].rank0(index);
      offset = rank_support_[level].rank0(offset);
    }
    return rank_matrix(level + 1, symbol, index, offset);
  }

  inline uint64_t select_tree(size_t const level, uint8_t const symbol,
    size_t occ, size_t const start, size_t const end) const {

    if (level == ws_.levels()) { return occ; }

    size_t const left = rank_support_[level].rank0(start);
    size_t const right = rank_support_[level].rank0(end);

    size_t const shift_for_bit = ws_.levels() - level - 1;
    if ((symbol >> shift_for_bit) & uint8_t(1)) {
      std::cout << "BIT IS ONE" << std::endl;
      occ = select_tree(level + 1, symbol, occ, start + right - left, end);
      std::cout << "RET " << select1_support_[level].select(start - left + occ) - start << std::endl;
      return select1_support_[level].select(start - left + occ) - start;
    } else {
      std::cout << "BIT IS ZERO" << std::endl;
      occ = select_tree(level + 1, symbol, occ, start, start + right -left);
      std::cout << "RET " << select0_support_[level].select(left + occ) - start << std::endl;
      return select0_support_[level].select(left + occ) - start;
    }
  }

  inline uint64_t select_matrix(size_t const level, uint8_t const symbol,
    size_t occ, size_t offset) const {

    if (level == ws_.levels()) { return occ + offset; }

    size_t const shift_for_bit = ws_.levels() - level - 1;
    if ((symbol >> shift_for_bit) & uint8_t(1)) {
      offset = ws_.zeros()[level] + rank_support_[level].rank1(offset);
      occ = select_matrix(level + 1, symbol, occ, offset);
      if (level == 0) {
        return select1_support_[level].select(occ - ws_.zeros()[level]);
      } else {
        return select1_support_[level].select(occ - ws_.zeros()[level]) + 1;
      }
    } else {
      offset = rank_support_[level].rank0(offset);
      occ = select_matrix(level + 1, symbol, occ, offset);
      if (level == 0) {
        return select0_support_[level].select(occ);
      } else { return select0_support_[level].select(occ) + 1; }
    }
  }

private:
  wavelet_structure& ws_;
  std::vector<rank_support_type> rank_support_;
  std::vector<select0_support_type> select0_support_;
  std::vector<select1_support_type> select1_support_;

}; // class query_support

/******************************************************************************/
