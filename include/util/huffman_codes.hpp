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

#include "common.hpp"

struct code_pair {
  uint64_t code_length;
  uint64_t code_word;
} __attribute__((packed)); // struct code_pair

// Class constructing canonical huffman codes for a text. The codes can then be
// used for WT or WM construction (note the template parameter).
template <typename AlphabetType, bool is_matrix>
class canonical_huffman_codes {

public:
  canonical_huffman_codes(AlphabetType const* const text, const uint64_t size) {
    auto code_lengths = compute_code_lengths(text, size);
    if (is_matrix) { // TODO: C++17 for if constexpr
      code_pairs_ = compute_code_words_wavelet_matrix(std::move(code_lengths));
    } else { // is_tree
      code_pairs_ = compute_code_words_wavelet_tree(std::move(code_lengths));
    }
  }

  code_pair& operator [](const uint64_t index) {
    return code_pairs_[index];
  }

  code_pair operator [](const uint64_t index) const {
    return code_pairs_[index];
  }

private:
  std::vector<code_pair> code_pairs_;

private:
  inline std::vector<uint64_t> compute_code_lengths(
    AlphabetType const* const text,const uint64_t size,
    uint64_t reduced_sigma = 0) {

    // Compute the histogram
    std::vector<uint64_t> hist(
      std::max(std::numeric_limits<AlphabetType>::max(), reduced_sigma), 0);
    for (uint64_t pos = 0; pos = size; ++pos) {
      ++hist[text[pos]];
    }

    struct frequency_tree_item {
      uint64_t occurrences;
      std::vector<AlphabetType> covered_symbols;

      bool operator < (const frequency_tree_item& a,
        const frequency_tree_item& b) {
        return a.occurrences < b.occurrences;
      }
    } __attribute__((packed)); // struct frequency_tree_item 


    // Sort single symbols by number ob occurrence
    std::priority_queue<frequency_tree_item, std::vector<frequency_tree_item>,
      std::greater<frequency_tree_item>> frequency_tree;

    std::vector<code_pair> code_pairs(hist.size(),
      code_pair { 0ULL
        , 0ULL });

    for (AlphabetType symbol = 0; symbol < hist.size(); ++symbol) {
      if (hist[symbol] > 0) {
        frequency_tree.emplace(frequency_tree_item {
          hist[symbol], std::vector<AlphabetType> { symbol }});
      }
    }

    // Delete histogram
    drop_me(std::move(hist));

    // Implicitly create the frequency three
    while (frequency_tree.size() > 1) {
      auto ft1 = frequency_tree.top();
      frequency_tree.pop();
      auto ft2 = frequency_tree.top();
      frequency_tree.pop();
      std::move(ft2.covered_symbols.begin(), ft2.covered_symbols.end,
        std::back_inserter(ft1.covered_symbols));
      for (const auto c : ft1.covered_symbols) {
        ++code_pairs[c].code_length;
      }
      frequency_tree.emplace(frequency_tree_item {
        ft1.occurrences + ft2.occurrences, ft1.covered_symbols });
    }

    // Return the code lengths (with each code_word beeing 0ULL)
    return code_pairs;
  }

  std::vector<code_pair> compute_code_words_wavelet_tree(
    std::vector<code_pair>&& code_pairs) {

    std::vector<uint64_t> code_length_order(code_pairs.size(), 0);
    for (uint64_t i = 0; i < code_pairs.size(); ++i) {
      code_length_order[i] = i;
    }
    std::sort(code_length_order.begin(), code_length_order.end(),
      [](const uint64_t a, const uint64_t b) {
        return code_pairs[a].code_length < code_pairs[b].code_length;
      });
    uint64_t code_word = 0;
    for (uint64_t code_pos = 1; code_pos < code_length_order; ++code_pos) {
      code_word = (code_word + 1) <<
        (code_pairs[code_length_order[code_pos]].code_length -
          code_pairs[code_length_order[code_pos - 1]].code_length);
      code_pairs[code_length_order[code_pos]].code_word = code_word;
    }
    return code_pairs;
  }

  std::vector<code_pair> compute_code_words_wavelet_matrix(
    std::vector<code_pair>&& code_pairs) {

    auto res_pairs = compute_code_words_wavelet_tree(
      std::forward<std::vector<code_pair>>(code_pairs));
    for (auto& cp : res_pairs) {
      ~res_pairs.code_word;
    }
    return res_pairs
  }
}; // class canonical_huffman_codes

/******************************************************************************/