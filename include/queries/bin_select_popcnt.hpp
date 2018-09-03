/*******************************************************************************
 * include/queries/bin_select_popcnt.hpp
 *
 * Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

/* Based on
 * @inproceedings{Zhou2013RankSelect,
 *    author    = {Dong Zhou and David G. Andersen and Michael Kaminsky},
 *    title     = {Space-Efficient, High-Performance Rank and Select Structures
 *                 on Uncompressed Bit Sequences},
 *    booktitle = {12th International Symposium on Experimental Algorithms
 *                 ({SEA})},
 *    series    = {LNCS},
 *    volume    = {7933},
 *    pages     = {151--163},
 *    publisher = {Springer},
 *    year      = {2013},
 * }
 */

#pragma once

#include <string> 
#include <vector>

#include "queries/bin_rank_popcnt.hpp"

template <uint8_t select_value>
class bin_select_popcnt {
  // static_assert(select_value > 1, std::to_string(select_value));

public:
  bin_select_popcnt(bin_rank_popcnt& rank_support)
  : rank_support_(rank_support), l0_(rank_support_.l0_),
    l12_(rank_support_.l12_)  {

    // TODO: Add sampling!
  }

  bin_select_popcnt() = delete;
  bin_select_popcnt(bin_select_popcnt const&) = delete;
  bin_select_popcnt& operator =(bin_select_popcnt const&) = delete;
  bin_select_popcnt(bin_select_popcnt&&) = default;
  bin_select_popcnt& operator =(bin_select_popcnt&&) = default;

  inline size_t select(size_t occ) const {
    if (occ == 0) { return 0; }

    if constexpr (select_value == 0) { return select0(occ); }
    else { return select1(occ); }
  }

private:
  inline size_t select0(size_t occ) const {
    size_t l0_pos = 0;
    while (l0_pos + 1 < l0_.size() &&
      ((l0_pos + 1) * rank_support_.l0_bit_size_ ) - l0_[l0_pos + 1] < occ) {
      ++l0_pos;
    }
    occ -= (l0_pos * rank_support_.l0_bit_size_ ) - l0_[l0_pos];

    size_t l1_pos = 0;
    while (l1_pos + 1 < l12_.size() &&
      ((l1_pos + 1) * rank_support_.l1_bit_size_) -
      rank_support_.get_l1_entry(l12_[l1_pos + 1]) < occ) {
      ++ l1_pos;
    }
    occ -= (l1_pos * rank_support_.l1_bit_size_) -
      rank_support_.get_l1_entry(l12_[l1_pos]);

    size_t l2_pos = 0;
    while (l2_pos < 3 &&
      occ > (rank_support_.l2_bit_size_ -
        rank_support_.get_l2_entry(l12_[l1_pos], l2_pos))) {
      occ -= (rank_support_.l2_bit_size_ -
        rank_support_.get_l2_entry(l12_[l1_pos], l2_pos++));
    }

    size_t last_pos = (rank_support_.l2_block_cover_ * l2_pos) +
      (rank_support_.l1_block_cover_ * l1_pos) +
      (rank_support_.l0_block_cover_ * l0_pos);

    size_t additional_words = 0;
    size_t popcnt = 0;
    while ((popcnt = __builtin_popcountll(
      ~(*(rank_support_.data_ + last_pos + additional_words)))) < occ) {
      ++additional_words;
      occ -= popcnt;
    }

    uint64_t final_word =
      ~(*(rank_support_.data_ + last_pos + additional_words));

    size_t last_bits = 0;
    uint64_t const mask = (1ULL << 63);
    while (occ > 0) {
      if (final_word & mask) { --occ; }
      ++last_bits;
      final_word <<= 1;
    }

    last_bits -= (last_bits > 0) ? 1 : 0;

    return (rank_support_.l2_bit_size_ * l2_pos) +
      (rank_support_.l1_bit_size_ * l1_pos) +
      (rank_support_.l0_bit_size_ * l0_pos) +
      (additional_words * 64) + last_bits;
  }

  inline size_t select1(size_t occ) const {
    size_t l0_pos = 0;
    while (l0_pos + 1 < l0_.size() && l0_[l0_pos + 1] < occ) { ++l0_pos; }
    occ -= l0_[l0_pos];

    size_t l1_pos = 0;
    while (l1_pos + 1 < l12_.size() &&
      rank_support_.get_l1_entry(l12_[l1_pos + 1]) < occ) {
      ++ l1_pos;
    }
    occ -= rank_support_.get_l1_entry(l12_[l1_pos]);

    size_t l2_pos = 0;
    while (l2_pos < 3 &&
      occ > rank_support_.get_l2_entry(l12_[l1_pos], l2_pos)) {
      occ -= rank_support_.get_l2_entry(l12_[l1_pos], l2_pos++);
    }

    size_t last_pos = (rank_support_.l2_block_cover_ * l2_pos) +
      (rank_support_.l1_block_cover_ * l1_pos) +
      (rank_support_.l0_block_cover_ * l0_pos);

    size_t additional_words = 0;
    size_t popcnt = 0;
    while ((popcnt = __builtin_popcountll(
      *(rank_support_.data_ + last_pos + additional_words))) < occ) {
      ++additional_words;
      occ -= popcnt;
    }

    uint64_t final_word = *(rank_support_.data_ + last_pos + additional_words);

    size_t last_bits = 0;
    uint64_t const mask = (1ULL << 63);
    while (occ > 0) {
      if (final_word & mask) { --occ; }
      ++last_bits;
      final_word <<= 1;
    }

    last_bits -= (last_bits > 0) ? 1 : 0;

    return (rank_support_.l2_bit_size_ * l2_pos) +
      (rank_support_.l1_bit_size_ * l1_pos) +
      (rank_support_.l0_bit_size_ * l0_pos) +
      (additional_words * 64) + last_bits;
  }  

private:
  // static constexpr size_t sample_rate = 8192;

  bin_rank_popcnt& rank_support_;
  std::vector<uint64_t>& l0_;
  std::vector<uint64_t>& l12_;

  // std::vector<uint32_t> samples_;
}; // class bin_select_popcnt

using bin_select0_popcnt = bin_select_popcnt<0>;
using bin_select1_popcnt = bin_select_popcnt<1>;

/******************************************************************************/
