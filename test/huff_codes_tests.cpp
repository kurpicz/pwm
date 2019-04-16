/*******************************************************************************
 * test/huffman_code_tests.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>

#include "util/util.hpp"

#include "huffman/huff_codes.hpp"
#include "huffman/huff_level_sizes_builder.hpp"
#include "util/alphabet_util.hpp"
#include "util/common.hpp"
#include "util/debug.hpp"
#include "util/file_util.hpp"

TEST(huffman_code_tests, wt_codes) {
  test::roundtrip_batch([](const std::string& s) {
    auto text = std::vector<uint8_t>(s.begin(), s.end());
    histogram<uint8_t> hist { text.data(), text.size() };
    level_sizes_builder<uint8_t> builder { std::move(hist) };
    canonical_huff_codes<uint8_t, false> chc(builder);
    std::vector<uint64_t> encoded_text;
    std::vector<uint64_t> code_lengths;
    for (const auto c : text) {
      const auto code = chc.encode_symbol(c);
      encoded_text.emplace_back(code.code_word());
      code_lengths.emplace_back(code.code_length());
    }
    ASSERT_EQ(encoded_text.size(), text.size());
    for (uint64_t i = 0; i < encoded_text.size(); ++i) {
      EXPECT_EQ(chc.decode_symbol(code_lengths[i], encoded_text[i]), text[i])
        << "pos: " << i << ", length: " << code_lengths[i] << ", symbol: "
        << encoded_text[i];
    }
  });
}

TEST(huffman_code_tests, wt_codes_reduced) {
  test::roundtrip_batch([](const std::string& s) {
    auto text = std::vector<uint8_t>(s.begin(), s.end());
    uint64_t max_char = reduce_alphabet(text);
    uint64_t reduced_sigma = levels_for_max_char(max_char);
    histogram<uint8_t> hist { text.data(), text.size(), reduced_sigma };
    level_sizes_builder<uint8_t> builder { std::move(hist) };
    canonical_huff_codes<uint8_t, false> chc(builder);
    std::vector<uint64_t> encoded_text;
    std::vector<uint64_t> code_lengths;
    for (const auto c : text) {
      const auto code = chc.encode_symbol(c);
      encoded_text.emplace_back(code.code_word());
      code_lengths.emplace_back(code.code_length());
    }
    ASSERT_EQ(encoded_text.size(), text.size());
    for (uint64_t i = 0; i < encoded_text.size(); ++i) {
      EXPECT_EQ(chc.decode_symbol(code_lengths[i], encoded_text[i]), text[i])
        << "pos: " << i << ", length: " << code_lengths[i] << ", symbol: "
        << encoded_text[i];
    }
  });
}

TEST(huffman_code_tests, wm_codes) {
  test::roundtrip_batch([](const std::string& s) {
    auto text = std::vector<uint8_t>(s.begin(), s.end());
    histogram<uint8_t> hist { text.data(), text.size() };
    level_sizes_builder<uint8_t> builder { std::move(hist) };
    canonical_huff_codes<uint8_t, true> chc(builder);
    std::vector<uint64_t> encoded_text;
    std::vector<uint64_t> code_lengths;
    for (const auto c : text) {
      const auto code = chc.encode_symbol(c);
      encoded_text.emplace_back(code.code_word());
      code_lengths.emplace_back(code.code_length());
    }
    ASSERT_EQ(encoded_text.size(), text.size());
    for (uint64_t i = 0; i < encoded_text.size(); ++i) {
      EXPECT_EQ(chc.decode_symbol(code_lengths[i], encoded_text[i]), text[i])
        << "pos: " << i << ", length: " << code_lengths[i] << ", symbol: "
        << encoded_text[i];
    }
  });
}

TEST(huffman_code_tests, wm_codes_reduced) {
  test::roundtrip_batch([](const std::string& s) {
    auto text = std::vector<uint8_t>(s.begin(), s.end());
    uint64_t max_char = reduce_alphabet(text);
    uint64_t reduced_sigma = levels_for_max_char(max_char);
    histogram<uint8_t> hist { text.data(), text.size(), reduced_sigma };
    level_sizes_builder<uint8_t> builder { std::move(hist) };
    canonical_huff_codes<uint8_t, true> chc(builder);
    std::vector<uint64_t> encoded_text;
    std::vector<uint64_t> code_lengths;
    for (const auto c : text) {
      const auto code = chc.encode_symbol(c);
      encoded_text.emplace_back(code.code_word());
      code_lengths.emplace_back(code.code_length());
    }
    ASSERT_EQ(encoded_text.size(), text.size());
    for (uint64_t i = 0; i < encoded_text.size(); ++i) {
      EXPECT_EQ(chc.decode_symbol(code_lengths[i], encoded_text[i]), text[i])
        << "pos: " << i << ", length: " << code_lengths[i] << ", symbol: "
        << encoded_text[i];
    }
  });
}

/******************************************************************************/
