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

class bin_rank_popcnt {

public:
  bin_rank_popcnt(uint64_t const * data, const size_t size)
  : l0_(size / upper_block_cover_, 0ULL), l12_(size / l12_block_cover_, 0ULL),
    data_(data), size_(size) {

    auto cur_pos = data;
    const auto end_pos = data + size;

    size_t l0_pos = 0;
    size_t l12_pos = 0;
    uint32_t l1_entry = 0UL;

    while (cur_pos != end_pos) {
      l12_[l12_pos] = set_l1_entry(l12_[l12_pos], l1_entry);
      for (size_t i = 0; i < 3; ++i) {
        uint32_t l2_entry = 0UL;
        for (size_t j = 0; j < basic_block_cover_; ++j) {
          l2_entry += __builtin_popcountll(*(cur_pos++));
        }
        l1_entry += l2_entry;
        l12_[l12_pos] = set_l2_entry(l12_[l12_pos], l2_entry, i);        
      }
      for (size_t j = 0; j < basic_block_cover_; ++j) {
        l1_entry += __builtin_popcountll(*(cur_pos++));
      }
      ++l12_pos;
      if (l12_pos % upper_block_cover_ == 0) {
        if (l0_pos > 0) { l0_[l0_pos] += l0_[l0_pos - 1]; }
        l0_[l0_pos++] += l1_entry;
        l1_entry = 0UL;
      }
    }
  }

  bin_rank_popcnt(bin_rank_popcnt const&) = delete;
  bin_rank_popcnt& operator =(bin_rank_popcnt const&) = delete;
  bin_rank_popcnt(bin_rank_popcnt&&) = default;
  bin_rank_popcnt& operator =(bin_rank_popcnt&&) = default;


  inline size_t rank1(const size_t index) const {
    size_t result = 0;
    size_t remaining_bits = index;
    // Find L0 block
    size_t l0_pos = remaining_bits / upper_block_bit_size_;
    remaining_bits -= (l0_pos * upper_block_bit_size_);
    if (l0_pos > 0) { result += l0_[l0_pos - 1]; }

    // Find L1/L2 block
    size_t l1_pos = remaining_bits / l12_block_bit_size_;
    remaining_bits -= (l1_pos * l12_block_bit_size_);
    const uint64_t l12_block = l12_[l1_pos];
    result += get_l1_entry(l12_block);

    size_t l2_pos = remaining_bits / basic_block_bit_size_;
    remaining_bits -= (l2_pos * basic_block_bit_size_);
    for (size_t i = 0; i < l2_pos; ++i) {
      result += get_l2_entry(l12_block, i);
    }

    uint64_t const * remaining_data = data_ + (l1_pos * l12_block_cover_) +
      (l2_pos * basic_block_cover_);
    while (remaining_bits >= 64) {
      result += __builtin_popcountll(*(remaining_data++));
      remaining_bits -= 64;
    }
    if (remaining_bits > 0) {
      result += __builtin_popcountll(*remaining_data >> (64 - remaining_bits));
    }
    return result;
  }

  inline size_t rank0(size_t index) const {
    return index - rank1(index);
  }

private:
  // Note that we can only set the value once. If we want to set the value
  // multiple times, we first have to set the first 32 bits to 0.
  inline uint64_t set_l1_entry(const uint64_t word,
    const uint32_t l1_entry) const {
    uint64_t l1_entry_cast = static_cast<uint64_t>(l1_entry);
    return word | (l1_entry_cast << 32); 
  }

  inline uint32_t get_l1_entry(const uint64_t word) const {
    return word >> 32;
  }

  // Note that we can only set each value once. If we want to set any value
  // multiple times, we firts ahve to the the corresponding bits to 0.
  inline uint64_t set_l2_entry(uint64_t word, const uint32_t l2_entry,
    const size_t entry_pos) const {
    uint64_t l2_entry_cast = static_cast<uint64_t>(l2_entry);
    return word | (l2_entry_cast << (22 - (entry_pos * 10)));
  }

  inline uint32_t get_l2_entry(const uint64_t word,
    const size_t entry_pos) const {
    return (word & ((~0ULL) >> (32 + (entry_pos * 10)))) >>
      (22 - (entry_pos * 10));
  }

private:
  static constexpr size_t bits_per_word_ = sizeof(uint64_t) * 8;
  static constexpr size_t basic_block_bit_size_ = 512;
  static constexpr size_t l12_block_bit_size_ = 4 * basic_block_bit_size_;
  static constexpr size_t upper_block_bit_size_ =
    std::numeric_limits<uint32_t>::max();

  static constexpr size_t basic_block_cover_ =
    basic_block_bit_size_ / bits_per_word_;
  static constexpr size_t l12_block_cover_ =
    l12_block_bit_size_ / bits_per_word_;
  static constexpr size_t upper_block_cover_ =
    upper_block_bit_size_ / bits_per_word_;

  std::vector<uint64_t> l0_;
  std::vector<uint64_t> l12_;
  uint64_t const * const data_;
  const size_t size_;
}; // class bin_rank_popcnt

/******************************************************************************/
