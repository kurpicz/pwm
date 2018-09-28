/*******************************************************************************
 * include/huffman/huff_codes.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <cstdint>
#include <map>
#include <queue>
#include <type_traits>

#include "util/common.hpp"
#include "util/macros.hpp"
#include "util/histogram.hpp"
#include "util/debug_assert.hpp"

struct code_pair {
  uint64_t code_length;
  uint64_t code_word;

  // Get the size-length prefix. This is necessary as shifting may not yield the
  // correct result due to the different length of the code words
  uint64_t prefix(const uint64_t size) const {
    DCHECK(size <= code_length);
    return code_word >> (code_length - size);
  }

  // Get the index-th bit of the code word. Note that we cannot simply shift
  // all the code words where they are used, as they have different lengths.
  bool operator [](const uint64_t index) const {
    DCHECK(index < code_length);
    return (code_word >> (code_length - index - 1)) & 1ULL;
  }

  bool operator ==(const code_pair& other) const {
    return std::tie(code_length, code_word) == std::tie(
      other.code_length, other.code_word);
  }

  bool operator <(const code_pair& other) const {
    return std::tie(code_length, code_word) < std::tie(
      other.code_length, other.code_word);
  }

  friend std::ostream& operator <<(std::ostream& os, const code_pair& cp) {
    os << "[ "
       << "length:" << cp.code_length
       << ", word: " << cp.code_word
       << ", bits: ";
    for(size_t i = 0; i < cp.code_length; i++) {
      os << int(cp[i]);
    }
    os << " ]";
    DCHECK((cp.code_word >> cp.code_length) == 0);
    return os;
  }
} PWM_ATTRIBUTE_PACKED; // struct code_pair

// Class constructing canonical huffman codes for a text. The codes can then be
// used for WT or WM construction (note the template parameter).
template <typename AlphabetType, bool is_tree>
class canonical_huff_codes {

public:
  canonical_huff_codes() = default;

  template<typename level_sizes_builder>
  canonical_huff_codes(level_sizes_builder& level_sizes)
  {
    compute_codes(level_sizes);
  }

  // Returns code_length and code_word for a given symbol, w.r.t. the text that
  // was used to create the canonical_huff_codes-instance.
  inline code_pair encode_symbol(const AlphabetType symbol) const {
    return code_pairs_[symbol];
  }

  uint64_t code_length(const AlphabetType symbol) const {
    return code_pairs_[symbol].code_length;
  }

  uint64_t code_word(const AlphabetType symbol) const {
    return code_pairs_[symbol].code_word;
  }

  // This for TESTING purposes ONLY!
  inline AlphabetType decode_symbol(const uint64_t length,
    const uint64_t encoded_symbol) const {

    return decode_table_.find({ length, encoded_symbol })->second;
  }

  inline const std::vector<code_pair>& code_pairs() const {
    return code_pairs_;
  }

private:
  std::vector<code_pair> code_pairs_;
  std::map<code_pair, AlphabetType> decode_table_;
private:

  template<typename level_sizes_builder>
  void compute_codes(level_sizes_builder& level_sizes)
  {
    const auto& histogram = level_sizes.get_histogram();

    struct frequency_tree_item {
      uint64_t occurrences;
      std::vector<AlphabetType> covered_symbols;

      bool operator > (const frequency_tree_item& other) const {
        return occurrences > other.occurrences;
      }
    }; // struct frequency_tree_item

    if (histogram.size() == 0) {
      return;
    }

    // Sort single symbols by number ob occurrence
    std::priority_queue<frequency_tree_item, std::vector<frequency_tree_item>,
      std::greater<frequency_tree_item>> frequency_tree;

    code_pairs_ = std::vector<code_pair>(histogram.max_symbol() + 1,
      code_pair { 0ULL, 0ULL });

    for (uint64_t i = 0; i < histogram.size(); ++i) {
      frequency_tree.emplace(frequency_tree_item {
        histogram[i].frequency,
        std::vector<AlphabetType> { histogram[i].symbol }});
    }

    // Corner case: Text consists of just one character
    if (PWM_UNLIKELY(frequency_tree.size() == 1)) {
      ++code_pairs_[frequency_tree.top().covered_symbols.front()].code_length;
    }

    // Implicitly create the frequency tree
    while (frequency_tree.size() > 1) {
      auto ft1 = frequency_tree.top();
      frequency_tree.pop();
      auto ft2 = frequency_tree.top();
      frequency_tree.pop();
      std::move(ft2.covered_symbols.begin(), ft2.covered_symbols.end(),
        std::back_inserter(ft1.covered_symbols));
      for (const auto c : ft1.covered_symbols) {
        ++code_pairs_[c].code_length;
      }
      frequency_tree.emplace(frequency_tree_item {
        ft1.occurrences + ft2.occurrences, ft1.covered_symbols });
    }

    std::vector<uint64_t> code_length_order(code_pairs_.size(), 0);
    for (uint64_t i = 0; i < code_pairs_.size(); ++i) {
      code_length_order[i] = i;
    }
    std::sort(code_length_order.begin(), code_length_order.end(),
      [&](const uint64_t a, const uint64_t b) {
        return code_pairs_[a].code_length < code_pairs_[b].code_length;
      });

    level_sizes.allocate_levels(
      code_pairs_[code_length_order.back()].code_length);

    if constexpr (!is_tree) {
      uint64_t code_nr = 0;
      uint64_t cur_length = 1;
      std::vector<uint64_t> code_words { 0ULL, 1ULL };

      while (code_nr < code_pairs_.size() &&
        code_pairs_[code_length_order[code_nr]].code_length == 0) { ++code_nr; }

      while (code_nr < code_pairs_.size()) {
        if (code_pairs_[code_length_order[code_nr]].code_length > cur_length) {
          for (uint i = cur_length;
            i < code_pairs_[code_length_order[code_nr]].code_length; ++i) {
            std::vector<uint64_t> new_code_words;

            new_code_words.reserve(code_words.size() << 1);
            for (const auto cw : code_words) {
              new_code_words.emplace_back(cw << 1);
            }
            for (const auto cw : code_words) {
              new_code_words.emplace_back((cw << 1) + 1);
            }
            code_words = std::move(new_code_words);
          }
          cur_length = code_pairs_[code_length_order[code_nr]].code_length;
        }
        while (code_nr < code_pairs_.size() &&
          code_pairs_[code_length_order[code_nr]].code_length == cur_length) {
          const uint64_t cur_code_pos = code_length_order[code_nr++];
          code_pairs_[cur_code_pos].code_word = code_words.back();
          code_words.pop_back();

          level_sizes.count(code_pairs_[cur_code_pos].code_length - 1,
                           cur_code_pos);
          decode_table_.emplace(std::make_pair(
            code_pairs_[cur_code_pos],
            AlphabetType(cur_code_pos)));
        }
      }
    } else { // if is_tree
      uint64_t code_word = 0ULL;
      uint64_t code_nr = 0;
      // The code lengths are correct, move to the second code word that has a
      // code_length > 0. The first one gets code_word = 0ULL.
      while (code_nr < code_pairs_.size() &&
        code_pairs_[code_length_order[code_nr]].code_length == 0) { ++code_nr; }

      auto generate_next_code = [&]{
        const uint64_t cur_code_pos = code_length_order[code_nr];
        auto& cur_code_pair = code_pairs_[cur_code_pos];

        // We use the bitwise negated word to ensure that a level in the WT is
        // cut of to the right. We also set all bits that do not belong to the
        // code word to 0. This helps when we decode naively when testing.

        cur_code_pair.code_word = (~code_word) &
          ((1ULL << cur_code_pair.code_length) - 1);
        decode_table_.emplace(std::make_pair(
          cur_code_pair, cur_code_pos));

        // Count the number of symbols that occur for each code length
        level_sizes.count(cur_code_pair.code_length - 1,
                         cur_code_pos);
      };

      generate_next_code();
      for (++code_nr; code_nr < code_pairs_.size(); ++code_nr) {
        // Create new code word
        code_word = (code_word + 1) << (
            code_pairs_[code_length_order[code_nr]].code_length -
            code_pairs_[code_length_order[code_nr - 1]].code_length);

        generate_next_code();
      }
    }
    level_sizes.drop_hist_and_finalize_level_sizes();
  }
}; // class canonical_huff_codes

/******************************************************************************/
