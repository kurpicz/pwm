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

static constexpr bool BRANCHLESS_SELECT = true;

public:
  bin_select_popcnt(bin_rank_popcnt& rank_support)
  : rank_support_(rank_support), l0_(rank_support_.l0_),
    l12_(rank_support_.l12_) {

    uint64_t const * const data = rank_support_.data_;
    size_t const size = rank_support_.size_;

    uint64_t cur_count = 0;
    uint64_t next_count = 0;
    uint64_t next_sample = 1;
    for (size_t i = 0; i < size; ++i) {
      uint64_t cur_value = data[i];

      if constexpr (select_value == 0) {
        cur_value = ~cur_value;
      }

      next_count = cur_count + __builtin_popcountll(cur_value);

      if (next_count > next_sample) {
        samples_.emplace_back((i * 64) + select1(cur_value, next_sample - cur_count));
        next_sample += SAMPLE_RATE_;
      }
      cur_count = next_count;
    }
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
    //*
    size_t l0_pos = 0;
    while (l0_pos + 1 < l0_.size() &&
      ((l0_pos + 1) * rank_support_.l0_bit_size_ ) - l0_[l0_pos + 1] < occ) {
      ++l0_pos;
    }
    occ -= (l0_pos * rank_support_.l0_bit_size_ ) - l0_[l0_pos];
    size_t l1_pos = 0;
    /*/

    size_t sample_pos = ((occ - 1) / SAMPLE_RATE_);
    
    if (occ % SAMPLE_RATE_ == 0) {
      return samples_[sample_pos];
    }
    occ %= SAMPLE_RATE_;

    size_t bit_pos = samples_.size() > 0 ? samples_[sample_pos] : 0;
    size_t l0_pos = bit_pos / rank_support_.l0_bit_size_;
    size_t l1_pos = bit_pos / rank_support_.l1_bit_size_;
    //*/
    while (l1_pos + 1 < l12_.size() &&
           ((l1_pos + 1) * rank_support_.l1_bit_size_) -
           l12_[l1_pos + 1].l1_value < occ) {
      ++ l1_pos;
    }
    occ -= (l1_pos * rank_support_.l1_bit_size_) - l12_[l1_pos].l1_value;

    size_t l2_pos = 0;
    while (l2_pos < 3 &&
      occ > (rank_support_.l2_bit_size_ - l12_[l1_pos][l2_pos])) {
      occ -= (rank_support_.l2_bit_size_ - l12_[l1_pos][l2_pos++]);
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
    while (l1_pos + 1 < l12_.size() && l12_[l1_pos + 1].l1_value < occ) {
      ++ l1_pos;
    }
    occ -= l12_[l1_pos].l1_value;

    size_t l2_pos = 0;
    while (l2_pos < 3 && occ > l12_[l1_pos][l2_pos]) {
      occ -= l12_[l1_pos][l2_pos++];
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

  inline uint32_t select1(uint64_t data, uint32_t rank) {
    // Select the bit position
    // From: https://graphics.stanford.edu/~seander/bithacks.html
    uint64_t v = data;     // Input value to find position with rank r
    unsigned int r = rank; // Input: bit's desired rank
    unsigned int s;        // Output: Resulting position of bit with rank r
    uint64_t a, b, c, d;   // Intermediate temporaries for bit count
    unsigned int t;        // Bit count temporary
    a =  v - ((v >> 1) & ~0UL/3);
    b = (a & ~0UL/5) + ((a >> 2) & ~0UL/5);
    c = (b + (b >> 4)) & ~0UL/0x11;
    d = (c + (c >> 8)) & ~0UL/0x101;
    t = (d >> 32) + (d >> 48);
    s = 64;
    if constexpr (BRANCHLESS_SELECT) {
      s -= ((t - r) & 256) >> 3; r -= (t & ((t - r) >> 8));
      t  = (d >> (s - 16)) & 0xff;
      s -= ((t - r) & 256) >> 4; r -= (t & ((t - r) >> 8));
      t  = (c >> (s - 8)) & 0xf;
      s -= ((t - r) & 256) >> 5; r -= (t & ((t - r) >> 8));
      t  = (b >> (s - 4)) & 0x7;
      s -= ((t - r) & 256) >> 6; r -= (t & ((t - r) >> 8));
      t  = (a >> (s - 2)) & 0x3;
      s -= ((t - r) & 256) >> 7; r -= (t & ((t - r) >> 8));
      t  = (v >> (s - 1)) & 0x1;
      s -= ((t - r) & 256) >> 8;
    } else {
      if (r > t) {s -= 32; r -= t;}
      t  = (d >> (s - 16)) & 0xff;
      if (r > t) {s -= 16; r -= t;}
      t  = (c >> (s - 8)) & 0xf;
      if (r > t) {s -= 8; r -= t;}
      t  = (b >> (s - 4)) & 0x7;
      if (r > t) {s -= 4; r -= t;}
      t  = (a >> (s - 2)) & 0x3;
      if (r > t) {s -= 2; r -= t;}
      t  = (v >> (s - 1)) & 0x1;
      if (r > t) s--;
    }
    // s = 65 - s;
    return 65 - s;
  }

private:
  static constexpr size_t SAMPLE_RATE_ = 8192;

  bin_rank_popcnt& rank_support_;
  std::vector<uint64_t>& l0_;
  std::vector<l12_entry>& l12_;

  std::vector<uint64_t> samples_;
}; // class bin_select_popcnt

using bin_select0_popcnt = bin_select_popcnt<0>;
using bin_select1_popcnt = bin_select_popcnt<1>;

/******************************************************************************/
