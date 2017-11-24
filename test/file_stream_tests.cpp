/*******************************************************************************
 * test/file_stream_tests.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>

#include <cstdio>

#include "test/util.hpp"

#include "util/file_stream.hpp"
#include "util/file_util.hpp"

TEST(ifile_stream, smoketest) {
  test::roundtrip_batch([&](const std::string& s){
    auto vec = std::vector<uint8_t>(s.begin(), s.end());
    std::vector<uint8_t> result(vec.size());
    std::string file_name = "test_file";
    vector_to_file(vec, file_name);

    {
      ifile_stream<uint8_t> fs(file_name);
      for (uint64_t i = 0; i < vec.size(); ++i) {
        result[i] = fs[i];
      }
    }
    if (remove(file_name.c_str()) != 0) {
      ASSERT_TRUE(false) << "Could not remove file, something went wrong";
    }
    for (uint64_t i = 0; i < vec.size(); ++i) {
      ASSERT_EQ(vec[i], result[i]) << "Failure when reading vector from file at"
                                   << " position " << i;
    }
  });
}

TEST(ofile_stream, smoketest) {
  test::roundtrip_batch([&](const std::string& s){
    auto vec = std::vector<uint8_t>(s.begin(), s.end());
    std::string file_name = "test_file";
    {
      ofile_stream<uint8_t> fs(file_name);
      for (uint64_t i = 0; i < vec.size(); ++i) {
        fs[i] = vec[i];
      }
    }
    auto result = file_to_vector<1>(file_name);

    if (remove(file_name.c_str()) != 0) {
      ASSERT_TRUE(false) << "Could not remove file, something went wrong";
    }
    for (uint64_t i = 0; i < vec.size(); ++i) {
      ASSERT_EQ(vec[i], result[i]) << "Failure when reading vector from file at"
                                   << " position " << i;
    }
  });
}

/******************************************************************************/
