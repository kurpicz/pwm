/*******************************************************************************
 * include/util/alphabet_util.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <limits>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "external_memory/stxxl_helper.hpp"

template <typename AlphabetType>
static uint64_t reduce_alphabet(std::vector<AlphabetType>& text) {
  uint64_t max_char = uint64_t(0);
  if constexpr (std::is_same<AlphabetType, uint8_t>::value) {
    std::array<uint64_t, std::numeric_limits<uint8_t>::max()> occ;
    occ.fill(0);
    for (auto& c : text) {
      if (occ[c] == 0) {
        occ[c] = ++max_char;
      }
      c = occ[c] - 1;
    }
    --max_char;
  } else {
    std::unordered_map<AlphabetType, uint64_t> word_list;
    for (const AlphabetType& character : text) {
      auto result = word_list.find(character);
      if (result == word_list.end()) {
        word_list.emplace(character, max_char++);
      }
    }
    --max_char;
    for (uint64_t i = 0; i < text.size(); ++i) {
      text[i] = static_cast<AlphabetType>(word_list.find(text[i])->second);
    }
  }
  return max_char;
}

template <typename AlphabetType>
static uint64_t reduce_alphabet(const stxxlvector<AlphabetType>& text,
                                stxxlvector<AlphabetType>& result,
                                uint64_t text_size) {
  if(text_size <= 0) text_size = text.size();
  else text_size = std::min(text_size, uint64_t(text.size()));

  uint64_t max_char = uint64_t(0);
  result.resize(0);
  result.reserve(text_size);
  stxxlreader<AlphabetType> reader(text);
  stxxlwriter<AlphabetType> writer(result);
  uint64_t counter = 0;

  if constexpr (std::is_same<AlphabetType, uint8_t>::value) {
    std::array<uint64_t, std::numeric_limits<uint8_t>::max()> occ;
    occ.fill(0);
    for (const auto& c : reader) {
      if (occ[c] == 0) {
        occ[c] = ++max_char;
      }
      writer << occ[c] - 1;
      if(++counter == text_size) break;
    }
    --max_char;
  } else {
    std::unordered_map<AlphabetType, uint64_t> word_list;
    for (const auto& c : reader) {
      auto result = word_list.find(c);
      if (result == word_list.end()) {
        word_list.emplace(c, max_char++);
      }
      writer << static_cast<AlphabetType>(word_list.find(c)->second);
      if(++counter == text_size) break;
    }
    --max_char;
  }
  writer.finish();
  return max_char;
}

[[maybe_unused]] static uint64_t levels_for_max_char(uint64_t max_char) {
  uint64_t levels = 0;
  while (max_char) {
    max_char >>= 1;
    ++levels;
  }
  return levels;
}

template <typename AlphabetType>
static uint64_t no_reduction_alphabet(const std::vector<AlphabetType>& /*t*/) {
  uint64_t max_char = std::numeric_limits<AlphabetType>::max();
  uint64_t levels = 0;
  while (max_char) {
    max_char >>= 1;
    ++levels;
  }
  return levels;
}

/******************************************************************************/
