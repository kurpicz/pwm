/*******************************************************************************
 * test/construction_test.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>

#include "util/util.hpp"

#include "benchmark/algorithm.hpp"
#include "huffman/huff_codes.hpp"
#include "huffman/huff_decode.hpp"
#include "huffman/huff_level_sizes_builder.hpp"
#include "util/alphabet_util.hpp"
#include "util/common.hpp"
#include "util/decode.hpp"
#include "util/file_util.hpp"

TEST(wavelet, no_alphabet_reduction) {
  auto& algo_list = algorithm_list::get_algorithm_list();
  for (const auto& a : algo_list) {
    if (a->word_width() == 1 && !a->is_huffman_shaped()) {
      a->print_info();
      test::roundtrip_batch([&](const std::string& s){
        auto vec = std::vector<uint8_t>(s.begin(), s.end());
        uint64_t levels = no_reduction_alphabet(vec);
        wavelet_structure bvz = a->compute_bitvector(&vec, vec.size() , levels);
        if (a->is_tree()) {
          auto decoded_s = decode_wt(bvz.bvs(), vec.size());
          ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
        } else {
          auto decoded_s = decode_wm(bvz.bvs(), bvz.zeros(), vec.size());
          ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
        }
      });
    }
  }
}

/*
TODO: Broken due to levels == 0 corner case

TEST(construction, wavelet_alphabet_reduction) {
  auto& algo_list = algorithm_list::get_algorithm_list();
  for (const auto& a : algo_list) {
    if (a->word_width() == 1 && !a->is_huffman_shaped()) {
      a->print_info();
      test::roundtrip_batch([&](const std::string& s){
        auto vec = std::vector<uint8_t>(s.begin(), s.end());
        auto max_char = reduce_alphabet(vec);
        uint64_t levels = levels_for_max_char(max_char);
        wavelet_structure bvz = a->compute_bitvector(&vec, vec.size() , levels);
        if (a->is_tree()) {
          auto decoded_s = decode_wt(bvz.bvs(), vec.size());
          auto decoded_vec = std::vector<uint8_t>(decoded_s.begin(), decoded_s.end());
          ASSERT_EQ(vec, decoded_vec) << "Failure (Algorithm): " << a->name();
        } else {
          auto decoded_s = decode_wm(bvz.bvs(), bvz.zeros(), vec.size());
          auto decoded_vec = std::vector<uint8_t>(decoded_s.begin(), decoded_s.end());
          ASSERT_EQ(vec, decoded_vec) << "Failure (Algorithm): " << a->name();
        }
      });
    }
  }
}
*/

TEST(huffman_shaped_wavelet, alphabet_reduction) {
  auto& algo_list = algorithm_list::get_algorithm_list();
  for (const auto& a : algo_list) {
    if (a->word_width() == 1 && a->is_huffman_shaped()) {
      a->print_info();
      test::roundtrip_batch([&](const std::string& s){
        auto vec = std::vector<uint8_t>(s.begin(), s.end());
        auto max_char = reduce_alphabet(vec);
        uint64_t levels = levels_for_max_char(max_char);
        auto bvz = a->compute_bitvector(&vec, vec.size() , levels);
        histogram<uint8_t> hist { vec.data(), vec.size() };
        level_sizes_builder<uint8_t> builder { std::move(hist) };
        if (a->is_tree()) {
          const auto codes =
            canonical_huff_codes<uint8_t, true>(builder);
          auto decoded_s = decode_wt_huff(bvz.bvs(), codes);
          auto decoded_vec = std::vector<uint8_t>(decoded_s.begin(), decoded_s.end());
          ASSERT_EQ(vec, decoded_vec) << "Failure (Algorithm): " << a->name();
        } else {
          const auto codes =
            canonical_huff_codes<uint8_t, false>(builder);
          auto decoded_s = decode_wm_huff(bvz.bvs(), codes);
          auto decoded_vec = std::vector<uint8_t>(decoded_s.begin(), decoded_s.end());
          ASSERT_EQ(vec, decoded_vec) << "Failure (Algorithm): " << a->name();
        }
      });
    }
  }
}

/******************************************************************************/
