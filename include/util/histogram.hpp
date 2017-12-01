/*******************************************************************************
 * include/util/histogram.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstdint>
#include <limits>
#include <type_traits>

template <typename AlphabetType>
struct histogram_entry {
  const AlphabetType symbol;
  uint64_t frequency;
}; // struct histogram_entry

template <typename AlphabetType>
class histogram {

public:
  histogram(const AlphabetType* text, const uint64_t size,
    const uint64_t reduced_sigma = 0) : max_symbol_(0){

    if (std::is_same<AlphabetType, uint8_t>::value) {
      const uint64_t max_char = std::max(
        static_cast<decltype(reduced_sigma)>(
          std::numeric_limits<uint8_t>::max() + 1), reduced_sigma);
      std::vector<uint64_t> hist(max_char, 0);
      for (uint64_t pos = 0; pos < size; ++pos) {
        const AlphabetType cur_char = text[pos];
        max_symbol_ = std::max(max_symbol_, cur_char);
        ++hist[cur_char];
      }
      for (uint64_t pos = 0; pos < hist.size(); ++pos) {
        if (hist[pos] > 0) {
          data_.emplace_back(
            histogram_entry<AlphabetType> { AlphabetType(pos), hist[pos] });
        }
      }
    } else {
      std::unordered_map<AlphabetType, uint64_t> symbol_list;
      for (uint64_t pos = 0; pos < size; ++pos) {
        const AlphabetType cur_char = text[pos];
        max_symbol_ = std::max(max_symbol_, cur_char);
        auto result = symbol_list.find(cur_char);
        if (result == symbol_list.end()) {
          symbol_list.emplace(cur_char, 1);
        } else {
          ++(result->second);
        }
      }
      for (const auto& symbol : symbol_list) {
        data_.emplace_back(
          histogram_entry<AlphabetType> { symbol.first, symbol.second });
      }
    }
  }

  AlphabetType max_symbol() const {
    return max_symbol_;
  };

  uint64_t size() const {
    return data_.size();
  }

  histogram_entry<AlphabetType>& operator [](const uint64_t index) {
    return data_[index];
  }

  histogram_entry<AlphabetType> operator [](const uint64_t index) const {
    return data_[index];
  }

private:
  AlphabetType max_symbol_;
  std::vector<histogram_entry<AlphabetType>> data_;
}; // class histogram

/******************************************************************************/
