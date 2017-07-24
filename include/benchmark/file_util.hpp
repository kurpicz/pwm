/*******************************************************************************
 * include/file_util.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef FILE_UTIL_HEADER
#define FILE_UTIL_HEADER
 
#include <fstream>
#include <limits>
#include <unordered_map>
#include <vector>

#include "util/type_for_bytes.hpp"

template <uint8_t BytesPerWord>
static std::vector<typename type_for_bytes<BytesPerWord>::type> file_to_vector(
  const std::string& file_name) {
  std::ifstream stream(file_name.c_str(), std::ios::in | std::ios::binary);
  stream.seekg(0, std::ios::end);
  uint64_t size = stream.tellg() / BytesPerWord;
  stream.seekg(0);
  std::vector<typename type_for_bytes<BytesPerWord>::type> result(size);
  stream.read(reinterpret_cast<char*>(result.data()), size);
  stream.close();
  return result;
}

template <typename AlphabetType>
static uint64_t reduce_alphabet(std::vector<AlphabetType>& text) {
  std::unordered_map<AlphabetType, uint64_t> word_list;
  uint64_t max_char = 0;
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

#endif // FILE_UTIL_HEADER

/******************************************************************************/
