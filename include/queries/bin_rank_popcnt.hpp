/*******************************************************************************
 * include/queries/bin_rank_popcnt.hpp
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

#include <limits>
#include <vector>

#include "arrays/span.hpp"

struct l12_entry {
  l12_entry() = default;

  l12_entry(uint32_t l1, std::array<uint32_t, 3> const& l2) : l1_value(l1),
    l2_values(uint32_t(0b1111111111) & l2[0]) {
    l2_values |= ((uint32_t(0b1111111111) & l2[1]) << 10);
    l2_values |= ((uint32_t(0b1111111111) & l2[2]) << 20);
  }

  inline uint32_t operator [](size_t index) const {
    return (l2_values >> 10 * index) & 0b1111111111;
  }

  uint32_t l1_value;
  uint32_t l2_values;
}; //struct l12_entry

class bin_rank_popcnt {
  template <uint8_t select_value>
  friend class bin_select_popcnt;

public:
  bin_rank_popcnt(span<uint64_t const> data)
  : l0_((data.size() / l0_block_cover_) + 1, 0ULL),
    l12_((data.size() / l1_block_cover_) + 1),
    data_(data) {

    size_t l0_pos = 0;
    size_t l12_pos = 0;
    uint32_t l1_entry = 0UL;

    for (size_t pos = 0; pos < data_.size();) {
      uint32_t new_l1_entry = l1_entry;
      std::array<uint32_t, 3> l2_entries = { 0, 0, 0 };
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < l2_block_cover_ && pos < data_.size(); ++j) {
          l2_entries[i] += __builtin_popcountll(data_[pos++]);
        }
        new_l1_entry += l2_entries[i];
      }
      l12_[l12_pos++] = l12_entry(l1_entry, l2_entries);
      l1_entry = new_l1_entry;
      for (size_t j = 0; j < l2_block_cover_ && pos < data_.size(); ++j) {
        l1_entry += __builtin_popcountll(data_[pos++]);
      }

      if (l12_pos % l0_block_cover_ == 0) {
        if (l0_pos > 0) { l0_[l0_pos] += l0_[l0_pos - 1]; }
        l0_[l0_pos++] += l1_entry;
        l1_entry = 0;
      }
    }
  }

  bin_rank_popcnt() = delete;
  bin_rank_popcnt(bin_rank_popcnt const&) = delete;
  bin_rank_popcnt& operator =(bin_rank_popcnt const&) = delete;
  bin_rank_popcnt(bin_rank_popcnt&&) = default;
  bin_rank_popcnt& operator =(bin_rank_popcnt&&) = default;

  inline size_t rank0(size_t index) const {
    return index - rank1(index);
  }

  inline size_t rank1(const size_t index) const {
    size_t result = 0;
    size_t remaining_bits = index;
    // Find L0 block
    size_t l0_pos = remaining_bits / l0_bit_size_;
    remaining_bits -= (l0_pos * l0_bit_size_);
    if (l0_pos > 0) { result += l0_[l0_pos - 1]; }

    // Find L1/L2 block
    size_t l1_pos = remaining_bits / l1_bit_size_;
    remaining_bits -= (l1_pos * l1_bit_size_);
    const l12_entry l12 = l12_[l1_pos];
    result += l12.l1_value;

    size_t l2_pos = remaining_bits / l2_bit_size_;
    remaining_bits -= (l2_pos * l2_bit_size_);
    for (size_t i = 0; i < l2_pos; ++i) {
      result += l12[i];
    }

    size_t offset = (l1_pos * l1_block_cover_) + (l2_pos * l2_block_cover_);
    while (remaining_bits >= 64) {
      result += __builtin_popcountll(data_[offset++]);
      remaining_bits -= 64;
    }
    if (remaining_bits > 0) {
      result += __builtin_popcountll(data_[offset] >> (64 - remaining_bits));
    }
    return result;
  }

private:
  static constexpr size_t bits_per_word_ = sizeof(uint64_t) * 8;
  static constexpr size_t l2_bit_size_ = 512;
  static constexpr size_t l1_bit_size_ = 4 * l2_bit_size_;
  static constexpr size_t l0_bit_size_ =
    std::numeric_limits<uint32_t>::max();

  static constexpr size_t l2_block_cover_ = l2_bit_size_ / bits_per_word_;
  static constexpr size_t l1_block_cover_ = l1_bit_size_ / bits_per_word_;
  static constexpr size_t l0_block_cover_ = l0_bit_size_ / bits_per_word_;

  std::vector<uint64_t> l0_;
  std::vector<l12_entry> l12_;
  
  span<uint64_t const> data_;
}; // class bin_rank_popcnt

/******************************************************************************/
