/*******************************************************************************
 * test/huffman_codes_test.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>
#include "test/util.hpp"

#include "util/alphabet_util.hpp"
#include "util/common.hpp"
#include "util/debug.hpp"
#include "util/file_util.hpp"
#include "util/huffman_codes.hpp"

TEST(huffman_code_computation, smoketest) {
  test::roundtrip_batch([](const std::string& s) {
    auto text = std::vector<uint8_t>(s.begin(), s.end())
    
  });
  auto& algo_list = algorithm_list::get_algorithm_list();
  for (const auto& a : algo_list) {
    if (a->word_width() == 1) {
      a->print_info();
      test::roundtrip_batch([&](const std::string& s){
        auto vec = std::vector<uint8_t>(s.begin(), s.end());
        uint64_t levels = no_reduction_alphabet(vec);
        auto bvz = a->compute_bitvector(&vec, vec.size() , levels);
        if (a->is_tree()) {
          auto decoded_s = decode_wt(bvz.raw_bvs(), vec.size());
          ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
        } else {
          auto decoded_s = decode_wm(
            bvz.raw_bvs(), bvz.raw_zeros(), vec.size());
          ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
        }
      });
    }
  }
}

/******************************************************************************/
