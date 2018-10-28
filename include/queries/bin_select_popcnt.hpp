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
#include "queries/popcnt_traits.hpp"

template <uint8_t SelectValue>
class bin_select_popcnt {

static constexpr bool BRANCHLESS_SELECT = true;

public:
  bin_select_popcnt(bin_rank_popcnt& rank_support) : l0_(rank_support.l0_),
    l12_(rank_support.l12_), data_(rank_support.data_) {

    samples_offset_.push_back(0);

    uint64_t next_sample = 1;
    size_t l0_pos = 0;
    for (size_t l1_pos = 0; l1_pos + 1 < l12_.size();) {
      
      // Each L1 block contains at most one sample
      // TODO: We can skip 4 L1 blocks, because only then there can be another
      //       sample position.
      if constexpr (SelectValue == 0) {
        if (((l1_pos + 1) * popcnt_traits::l1_bit_size) -
            l12_[l1_pos + 1].l1_value > next_sample) {

          uint32_t remaining = next_sample -
            ((l1_pos * popcnt_traits::l1_bit_size) - l12_[l1_pos].l1_value);
          size_t l2_pos = 0;
          while (l2_pos < 3 && (popcnt_traits::l2_bit_size -
                                l12_[l1_pos][l2_pos]) < remaining) {
            remaining -= (popcnt_traits::l2_bit_size - l12_[l1_pos][l2_pos++]);
          }

          size_t word_pos = (l0_pos * popcnt_traits::l0_block_cover +
                             l1_pos * popcnt_traits::l1_block_cover +
                             l2_pos * popcnt_traits::l2_block_cover);

          while (word_pos < data_.size() &&
                 remaining > __builtin_popcountll(~data_[word_pos])) {
            remaining -= __builtin_popcountll(~data_[word_pos++]);
          }


          if (word_pos < data_.size()) {
            samples_.emplace_back(word_select1(~data_[word_pos], remaining) +
              64 * (word_pos - (l0_pos * popcnt_traits::l0_bit_size)));
            next_sample += SAMPLE_RATE_;
          }
        }
      } else {
        if (l12_[l1_pos + 1].l1_value > next_sample) {
          uint32_t remaining = next_sample - l12_[l1_pos].l1_value;
          size_t l2_pos = 0;
          while (l2_pos < 3 && l12_[l1_pos][l2_pos] < remaining) {
            remaining -= l12_[l1_pos][l2_pos++];
          }

          size_t word_pos = (l0_pos * popcnt_traits::l0_block_cover +
                             l1_pos * popcnt_traits::l1_block_cover +
                             l2_pos * popcnt_traits::l2_block_cover);

          while (word_pos < data_.size() &&
                 remaining > __builtin_popcountll(data_[word_pos])) {
            remaining -= __builtin_popcountll(data_[word_pos++]);
          }
          if (word_pos < data_.size()) {
            samples_.emplace_back(word_select1(data_[word_pos], remaining) +
              64 * (word_pos - (l0_pos * popcnt_traits::l0_bit_size)));
            next_sample += SAMPLE_RATE_;
          }
        }
      }

      if (++l1_pos % popcnt_traits::l0_block_cover == 0) {
        ++l0_pos;
        samples_offset_.push_back(samples_.size());
      }
    }
  }

  bin_select_popcnt() = delete;
  bin_select_popcnt(bin_select_popcnt const&) = delete;
  bin_select_popcnt& operator =(bin_select_popcnt const&) = delete;
  bin_select_popcnt(bin_select_popcnt&&) = default;
  bin_select_popcnt& operator =(bin_select_popcnt&&) = default;

  inline size_t select(size_t occ) const {
    if (occ == 0) { return 0; } 

    if constexpr (SelectValue == 0) { return select0(occ); }
    else { return select1(occ); }
  }

private:
  inline uint32_t word_select1(uint64_t data, uint32_t rank) const {
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
    return 64 - s; // We want the rank to be in [0..64)
  }

  inline size_t select0(size_t occ) const {
    size_t l0_pos = 0;
    while (l0_pos + 1 < l0_.size() && l0_[l0_pos + 1] < occ) { ++l0_pos; }
    occ -= l0_[l0_pos];

    size_t l1_pos = samples_[occ / SAMPLE_RATE_];
    if ((occ - 1) % SAMPLE_RATE_ == 0) { return l1_pos; }
    else { l1_pos /= popcnt_traits::l1_bit_size; }
    while (l1_pos + 1 < l12_.size() &&
           (((l1_pos + 1) * popcnt_traits::l1_bit_size) - l12_[l1_pos + 1].l1_value) < occ) {
      ++l1_pos;
    }

    occ -= (l1_pos * popcnt_traits::l1_bit_size) - l12_[l1_pos].l1_value;

    size_t l2_pos = 0;
    while (l2_pos < 3 && occ > popcnt_traits::l2_bit_size - l12_[l1_pos][l2_pos]) {
      occ -= popcnt_traits::l2_bit_size - l12_[l1_pos][l2_pos++];
    }

    size_t last_pos = (popcnt_traits::l2_block_cover * l2_pos) +
                      (popcnt_traits::l1_block_cover * l1_pos) +
                      (popcnt_traits::l0_block_cover * l0_pos);

    size_t additional_words = 0;
    size_t popcnt = 0;
    while ((popcnt =
           __builtin_popcountll(~data_[last_pos + additional_words])) < occ) {
      ++additional_words;
      occ -= popcnt;
    }

    return (popcnt_traits::l2_bit_size * l2_pos) +
      (popcnt_traits::l1_bit_size * l1_pos) +
      (popcnt_traits::l0_bit_size * l0_pos) +
      (additional_words * 64) +
      word_select1(~data_[last_pos + additional_words], uint32_t(occ)); 
  }

  inline size_t select1(size_t occ) const {
    size_t l0_pos = 0;
    while (l0_pos + 1 < l0_.size() && l0_[l0_pos + 1] < occ) { ++l0_pos; }
    occ -= l0_[l0_pos];

    size_t l1_pos = samples_[occ / SAMPLE_RATE_];
    if ((occ - 1) % SAMPLE_RATE_ == 0) { return l1_pos; }
    else { l1_pos /= popcnt_traits::l1_bit_size; }
    while (l1_pos + 1 < l12_.size() && l12_[l1_pos + 1].l1_value < occ) {
      ++l1_pos;
    }

    occ -= l12_[l1_pos].l1_value;

    size_t l2_pos = 0;
    while (l2_pos < 3 && occ > l12_[l1_pos][l2_pos]) {
      occ -= l12_[l1_pos][l2_pos++];
    }

    size_t last_pos = (popcnt_traits::l2_block_cover * l2_pos) +
                      (popcnt_traits::l1_block_cover * l1_pos) +
                      (popcnt_traits::l0_block_cover * l0_pos);

    size_t additional_words = 0;
    size_t popcnt = 0;
    while ((popcnt =
           __builtin_popcountll(data_[last_pos + additional_words])) < occ) {
      ++additional_words;
      occ -= popcnt;
    }

    return (popcnt_traits::l2_bit_size * l2_pos) +
      (popcnt_traits::l1_bit_size * l1_pos) +
      (popcnt_traits::l0_bit_size * l0_pos) +
      (additional_words * 64) +
      word_select1(data_[last_pos + additional_words], uint32_t(occ));
  }

private:
  static constexpr size_t SAMPLE_RATE_ = 8192;

  std::vector<uint64_t>& l0_;
  std::vector<l12_entry>& l12_;
  span<uint64_t const>& data_;

  std::vector<uint32_t> samples_;
  std::vector<uint64_t> samples_offset_;
}; // class bin_select_popcnt

using bin_select0_popcnt = bin_select_popcnt<0>;
using bin_select1_popcnt = bin_select_popcnt<1>;

/******************************************************************************/
