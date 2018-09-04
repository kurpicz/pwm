/*******************************************************************************
 * test/query_support_tests.cpp
 *
 * Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>
#include <vector>

#include "util/util.hpp"
#include "util/alphabet_util.hpp"
#include "util/print.hpp"

#include "queries/query_support.hpp"

#include "wx_naive.hpp"

TEST(access_tests, access_wt) {
  test::roundtrip_batch([&](std::string const& s){
    auto vec = std::vector<uint8_t>(s.begin(), s.end());
    uint64_t levels = no_reduction_alphabet(vec);

    auto wt = wx_naive<uint8_t, true>::compute(vec.data(),
      vec.size(), levels);

    query_support qs(wt);

    for (size_t i = 0; i < vec.size(); ++i) { ASSERT_EQ(qs.access(i), vec[i]); }
  });
}

TEST(access_tests, access_wm) {
  test::roundtrip_batch([&](std::string const& s){
    auto vec = std::vector<uint8_t>(s.begin(), s.end());
    uint64_t levels = no_reduction_alphabet(vec);

    auto wm = wx_naive<uint8_t, false>::compute(vec.data(),
      vec.size(), levels);

    query_support qs(wm);
    for (size_t i = 0; i < vec.size(); ++i) { ASSERT_EQ(qs.access(i), vec[i]); }
  });
}

TEST(rank_tests, rank_wt) {
  test::roundtrip_batch([&](std::string const& s){
    auto vec = std::vector<uint8_t>(s.begin(), s.end());
    uint64_t levels = no_reduction_alphabet(vec);

    auto wm = wx_naive<uint8_t, true>::compute(vec.data(),
      vec.size(), levels);

    std::vector<size_t> symbol_counts(256, 0);

    query_support qs(wm);
    for (size_t i = 0; i < vec.size(); ++i) {
      for (size_t j = 0; j < 256; ++j) {
        ASSERT_EQ(qs.rank(j, i), symbol_counts[j]);
      }
      ++symbol_counts[vec[i]];
    }
  });  
}

TEST(rank_tests, rank_wm) {
  test::roundtrip_batch([&](std::string const& s){
    auto vec = std::vector<uint8_t>(s.begin(), s.end());
    uint64_t levels = no_reduction_alphabet(vec);

    auto wm = wx_naive<uint8_t, false>::compute(vec.data(),
      vec.size(), levels);

    std::vector<size_t> symbol_counts(256, 0);

    query_support qs(wm);
    for (size_t i = 0; i < vec.size(); ++i) {
      for (size_t j = 0; j < 256; ++j) {
        ASSERT_EQ(qs.rank(j, i), symbol_counts[j]);
      }
      ++symbol_counts[vec[i]];
    }
  });
}

TEST(select_tests, select_wt) {
  test::roundtrip_batch([&](std::string const& s){
    auto vec = std::vector<uint8_t>(s.begin(), s.end());
    uint64_t levels = no_reduction_alphabet(vec);

    auto wt = wx_naive<uint8_t, true>::compute(vec.data(),
      vec.size(), levels);

    std::vector<size_t> symbol_counts(256, 0);

    query_support qs(wt);
    for (size_t i = 0; i < vec.size(); ++i) {
      EXPECT_EQ(qs.select(vec[i], ++symbol_counts[vec[i]]), i);
    }
  });
}

TEST(select_tests, select_wm) {
  test::roundtrip_batch([&](std::string const& s){
    auto vec = std::vector<uint8_t>(s.begin(), s.end());
    uint64_t levels = no_reduction_alphabet(vec);

    auto wm = wx_naive<uint8_t, false>::compute(vec.data(),
      vec.size(), levels);

    std::vector<size_t> symbol_counts(256, 0);

    query_support qs(wm);
    for (size_t i = 0; i < vec.size(); ++i) {
      ASSERT_EQ(qs.select(vec[i], ++symbol_counts[vec[i]]), i);
    }
  });
}

/******************************************************************************/
