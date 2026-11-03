/*******************************************************************************
 * test/util_tests.cpp
 *
 * Copyright (C) 2026 Donato Meoli
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "util/alphabet_util.hpp"
#include "util/file_util.hpp"

TEST(util_tests, reduce_alphabet_all_bytes) {
  // every byte value occurs, 255 first, so the reduced text is 0, 1, ..., 255
  std::vector<uint8_t> text(256);
  for (uint64_t i = 0; i < text.size(); ++i) {
    text[i] = static_cast<uint8_t>(255 - i);
  }
  text.push_back(255);
  text.push_back(0);
  const uint64_t max_char = reduce_alphabet(text);
  EXPECT_EQ(max_char, 255U);
  for (uint64_t i = 0; i < 256; ++i) {
    EXPECT_EQ(text[i], i);
  }
  EXPECT_EQ(text[256], 0U);
  EXPECT_EQ(text[257], 255U);
}

template <uint8_t BytesPerWord>
static void file_to_vector_roundtrip() {
  using word = typename type_for_bytes<BytesPerWord>::type;
  std::vector<word> words(1000);
  std::iota(words.begin(), words.end(), word(1));
  const std::string file_name =
      "util_tests_" + std::to_string(BytesPerWord) + ".bin";
  {
    std::ofstream out(file_name, std::ios::binary);
    out.write(reinterpret_cast<const char*>(words.data()),
              words.size() * sizeof(word));
  }
  EXPECT_EQ(file_to_vector<BytesPerWord>(file_name), words);
  const auto prefix = file_to_vector<BytesPerWord>(file_name, 10);
  EXPECT_EQ(prefix, std::vector<word>(words.begin(), words.begin() + 10));
  std::remove(file_name.c_str());
}

TEST(util_tests, file_to_vector_words) {
  file_to_vector_roundtrip<1>();
  file_to_vector_roundtrip<2>();
  file_to_vector_roundtrip<4>();
}

/******************************************************************************/
