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
    auto text = std::vector<uint8_t>(s.begin(), s.end());
    
  });
}

/******************************************************************************/
