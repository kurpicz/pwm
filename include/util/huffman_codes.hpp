/*******************************************************************************
 * include/util/huffman_codes.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <queue>
#include <type_traits>

#include "common.hpp"
#include "macros.hpp"

struct code_pair {
  uint64_t code_length;
  uint64_t code_word;
} PWM_ATTRIBUTE_PACKED; // struct code_pair

// Class constructing canonical huffman codes for a text. The codes can then be
// used for WT or WM construction (note the template parameter).
template <typename AlphabetType, bool is_matrix>
class canonical_huffman_codes {

public:
  canonical_huffman_codes(AlphabetType const* const text, const uint64_t size) {
    compute_code_lengths(text, size);
    compute_code_words();
    // if (is_matrix) { // TODO: C++17 for if constexpr
    //   static_assert(false,
    //     "NOT YET IMPLEMENTED. WE FIRST IMPLEMENTED IT FOR WTs");
    // }
  }

  // Returns code_length and code_word for a given symbol, w.r.t. the text that
  // was used to create the canonical_huffman_codes-instance.
  code_pair encode_symbol(AlphabetType symbol) const {
    return code_pairs_[symbol];
  }

  AlphabetType decode_symbol(const uint64_t encoded_symbol) const {
    return decode_table_[encoded_symbol];
  }

private:
  std::vector<code_pair> code_pairs_;
  std::vector<AlphabetType> decode_table_;

private:
  void compute_code_lengths(
    AlphabetType const* const text,const uint64_t size,
    uint64_t reduced_sigma = 0) {

    // Compute the histogram
    const uint64_t max_char = std::max(
      static_cast<decltype(reduced_sigma)>(
        std::numeric_limits<AlphabetType>::max()), reduced_sigma);
    std::vector<uint64_t> hist(max_char, 0);
    decode_table_ = std::vector<AlphabetType>(max_char);
    for (uint64_t pos = 0; pos < size; ++pos) {
      ++hist[text[pos]];
    }

    struct frequency_tree_item {
      uint64_t occurrences;
      std::vector<AlphabetType> covered_symbols;

      bool operator < (const frequency_tree_item& other) const {
        return occurrences < other.occurrences;
      }

      bool operator > (const frequency_tree_item& other) const {
        return occurrences > other.occurrences;
      }
    }; // struct frequency_tree_item 

    // Sort single symbols by number ob occurrence
    std::priority_queue<frequency_tree_item, std::vector<frequency_tree_item>,
      std::greater<frequency_tree_item>> frequency_tree;

    code_pairs_ = std::vector<code_pair>(hist.size(),
      code_pair { 0ULL, 0ULL });

    for (AlphabetType symbol = 0; symbol < hist.size(); ++symbol) {
      if (hist[symbol] > 0) {
        frequency_tree.emplace(frequency_tree_item {
          hist[symbol], std::vector<AlphabetType> { symbol }});
      }
    }

    // Cornder case: Text consists of just one character
    if (frequency_tree.size() == 1) {
      ++code_pairs_[frequency_tree.top().covered_symbols.front()].code_length;
    }

    // Implicitly create the frequency three
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
  }

  void compute_code_words() {
    std::vector<uint64_t> code_length_order(code_pairs_.size(), 0);
    for (uint64_t i = 0; i < code_pairs_.size(); ++i) {
      code_length_order[i] = i;
    }
    std::sort(code_length_order.begin(), code_length_order.end(),
      [&](const uint64_t a, const uint64_t b) {
        return code_pairs_[a].code_length < code_pairs_[b].code_length;
      });


    uint64_t code_word = 0;
    uint64_t code_nr = 0;
    // The code lengths are correct, move to the second code word that has a
    // code_lenght > 0. The first one gehts code_word = 0ULL.
    while (code_nr < code_pairs_.size() &&
      code_pairs_[code_length_order[code_nr++]].code_length == 0) { }
    decode_table_[0] = AlphabetType(code_length_order[code_nr - 1]);
    for (; code_nr < code_pairs_.size(); ++code_nr) {
      const uint64_t cur_code_pos = code_length_order[code_nr];

      // Create new code word
      code_word = (code_word + 1) << (code_pairs_[cur_code_pos].code_length -
              code_pairs_[code_length_order[code_nr - 1]].code_length);
      code_pairs_[cur_code_pos].code_word = code_word;
      decode_table_[code_word] = AlphabetType(cur_code_pos);
    }
  }
}; // class canonical_huffman_codes

/******************************************************************************/
