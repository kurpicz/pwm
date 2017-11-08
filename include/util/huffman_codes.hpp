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

template <typename AlphabetType>
struct decode_result {
  AlphabetType symbol;
  uint64_t code_length;
} __attribute__((packed));

// Class constructing canonical huffman codes for a text. The codes can then be
// used for WT or WM construction (note the template parameter).
template <typename AlphabetType, bool is_matrix>
class canonical_huffman_codes {

public:
  canonical_huffman_codes(AlphabetType const* const text, const uint64_t size) {
    compute_code_lengths(text, size);
    if (is_matrix) { // TODO: C++17 for if constexpr
      compute_code_words_wavelet_matrix();
    } else { // is_tree
      compute_code_words_wavelet_tree();
    }
  }

  // Returns code_length and code_word for a given symbol, w.r.t. the text that
  // was used to create the canonical_huffman_codes-instance.
  code_pair encode_symbol(AlphabetType symbol) const {
    return code_pairs_[symbol];
  }

  std::vector<uint64_t> encode_text(AlphabetType const* const text,
    const uint64_t size) {

    std::vector<uint64_t> result;

    uint64_t cur_word = 0ULL;
    uint64_t cur_word_pos = 0;
    for (uint64_t pos = 0, cur_word_pos = 0; pos < size; ++pos) {
      const auto code = encode_symbol(text[pos]);
      if (cur_word_pos + code.code_length > 64) { // split code_word
        cur_word << (64 - cur_word_pos);
        cur_word |= code.code_word >> (64 - cur_word_pos);
        result.emplace_back(cur_word);
        cur_word = 0ULL;
        cur_word |= (((1ULL << (64 - cur_word_pos)) - 1) | code.code_word);
        cur_word_pos = (64 - cur_word_pos);
      } else {
        cur_word << code.code_length;
        cur_word |= code.code_word;
        cur_word_pos += code.code_length;
      }
      result.emplace_back(cur_word)
    }
    return result;
  }

  // Given an array of 64-bit words and an offset (from the beginning of the
  // array), this function returns the first symbol (of the original alphabet)
  // and the length of the code_word that was used to encode the symbol.
  decode_result<AlphabetType> decode_symbol(
    uint64_t const* const encoded_vector, const uint64_t offset) {


  }

  std::vector<AlphabetType> decode_text(const std::vector<uint64_t>& enc_text) {

  }

private:
  std::vector<code_pair> code_pairs_;

private:
  void compute_code_lengths(
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

    code_pairs_ = std::vector<code_pair>(hist.size(),
      code_pair { 0ULL, 0ULL });

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
        ++code_pairs_[c].code_length;
      }
      frequency_tree.emplace(frequency_tree_item {
        ft1.occurrences + ft2.occurrences, ft1.covered_symbols });
    }
  }

  std::vector<code_pair> compute_code_words_wavelet_tree() {

    std::vector<uint64_t> code_length_order(code_pairs_.size(), 0);
    for (uint64_t i = 0; i < code_pairs_.size(); ++i) {
      code_length_order[i] = i;
    }
    std::sort(code_length_order.begin(), code_length_order.end(),
      [](const uint64_t a, const uint64_t b) {
        return code_pairs_[a].code_length < code_pairs_[b].code_length;
      });
    uint64_t code_word = 0;
    for (uint64_t code_pos = 1; code_pos < code_length_order; ++code_pos) {
      code_word = (code_word + 1) <<
        (code_pairs_[code_length_order[code_pos]].code_length -
          code_pairs_[code_length_order[code_pos - 1]].code_length);
      code_pairs_[code_length_order[code_pos]].code_word = code_word;
    }
  }

  std::vector<code_pair> compute_code_words_wavelet_matrix(
    std::vector<code_pair>&& code_pairs) {

    compute_code_words_wavelet_tree();
    for (auto& cp : code_pairs_) {
      ~code_pairs_.code_word;
    }
  }
}; // class canonical_huffman_codes

/******************************************************************************/
