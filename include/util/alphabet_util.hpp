/*******************************************************************************
 * include/util/alphabet_util.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <limits>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "util/file_stream.hpp"
#include "util/stxxl_helper.hpp"

template <typename AlphabetType>
static uint64_t reduce_alphabet(std::vector<AlphabetType>& text) {
  uint64_t max_char = uint64_t(0);
  if (std::is_same<AlphabetType, uint8_t>::value) { // TODO: C++17 if constexpr
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
    for (const auto& character : text) {
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
static uint64_t reduce_alphabet(const stxxlvector<AlphabetType>& text, stxxlvector<AlphabetType>& result) {
  uint64_t max_char = uint64_t(0);
  result.resize(0);
  result.reserve(text.size());
  stxxlreader<AlphabetType> reader(text);
  stxxlwriter<AlphabetType> writer(result);
  
  if (std::is_same<AlphabetType, uint8_t>::value) { // TODO: C++17 if constexpr
    std::array<uint64_t, std::numeric_limits<uint8_t>::max()> occ;
    occ.fill(0);
    for(const auto &c : reader) {
      if (occ[c] == 0) {
        occ[c] = ++max_char;
      }
      writer << occ[c] - 1;
    }
    --max_char;
  } else {
    std::unordered_map<AlphabetType, uint64_t> word_list;
    for(const auto &c : reader) {
      auto result = word_list.find(c);
      if (result == word_list.end()) {
        word_list.emplace(c, max_char++);
      }
      writer << static_cast<AlphabetType>(word_list.find(c)->second);
    }
    --max_char;
  }
  writer.finish();
  return max_char;
}

template <typename AlphabetType>
auto reduce_alphabet_stream(const std::string& file_name) {
  struct se_reduced_info {
    std::string file_name;
    uint64_t file_size;
    uint64_t max_char;
  }; // struct se_reduced_info

  ifile_stream<AlphabetType> ifs(file_name);
  const uint64_t file_size = ifs.file_size();
  const std::string reduced_file_name = file_name + std::string("_reduced");
  ofile_stream<AlphabetType> ofs(reduced_file_name);

  uint64_t max_char = 0;
  if (std::is_same<AlphabetType, uint8_t>::value) { // TODO: C++17 if constexpr
    std::array<uint64_t, std::numeric_limits<uint8_t>::max()> occ;
    occ.fill(0);
    for (uint64_t i = 0; i < file_size; ++i) {
      const AlphabetType c = ifs[i];
      if (occ[c] == 0) {
        occ[c] = ++max_char;
      }
      ofs[i] = occ[c] - 1;
    }
    --max_char;
  } else {
    std::unordered_map<AlphabetType, uint64_t> word_list;
    for (uint64_t i = 0; i < file_size; ++i) {
      const AlphabetType c = ifs[i];
      auto result = word_list.find(c);
      if (result == word_list.end()) {
        word_list.emplace(c, max_char++);
      }
    }
    ifs.reset_stream();
    --max_char;
    for (uint64_t i = 0; i < file_size; ++i) {
      ofs[i] = static_cast<AlphabetType>(word_list.find(ifs[i])->second);
    }
  }
  return  se_reduced_info { reduced_file_name, file_size, max_char };
}

[[gnu::unused]] // TODO: C++17 [[maybe_unused]] 
static uint64_t levels_for_max_char(uint64_t max_char) {
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
  return max_char;
}

/******************************************************************************/
