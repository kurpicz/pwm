/*******************************************************************************
 * test/query_support_tests.cpp
 *
 * Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>
#include <vector>

#include "test/util.hpp"
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

/******************************************************************************/
