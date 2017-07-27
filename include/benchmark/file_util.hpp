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

template<typename AlphabetType>
class Histogram  {
    std::unordered_map<AlphabetType, uint64_t> m_hist;
public:
    inline Histogram() = default;

    inline uint64_t& entry(AlphabetType chr) {
        return m_hist[chr];
    }
};

template<>
class Histogram<uint8_t>  {
    std::array<uint64_t, 256> m_hist;
public:
    inline Histogram() {
        m_hist.fill(0);
    }

    inline uint64_t& entry(uint8_t chr) {
        return m_hist[chr];
    }
};

template <typename AlphabetType>
static uint64_t reduce_alphabet(std::vector<AlphabetType>& text) {
  Histogram<AlphabetType> word_list;
  uint64_t max_char = 0;
  for (const AlphabetType& character : text) {
    auto& entry = word_list.entry(character);
    if (entry == 0) {
      entry = 1 + max_char++;
    }
  }
  --max_char;
  for (uint64_t i = 0; i < text.size(); ++i) {
    text[i] = static_cast<AlphabetType>(word_list.entry(text[i]) - 1);
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
