/*******************************************************************************
 * test/construction_test.cpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>

#include <random>
#include <bitset>

#include "util/util.hpp"
#include "util/debug.hpp"
#include "construction/pc_dd_fe/merge_external.hpp"

TEST(external_merge, random_bvs) {
  using vec_type = stxxlvector<uint64_t>;
  using reader_type = typename vec_type::bufreader_type;
  using writer_type = typename vec_type::bufwriter_type;

  const uint64_t test_count = 250;
  const uint64_t max_len = 250;

  std::random_device rd;
  std::mt19937_64 gen(rd());
  std::uniform_int_distribution<uint64_t> rnd_word;
  std::uniform_int_distribution<> rnd_length(0, max_len);

  for(uint64_t i = 0; i < test_count; ++i) {
    uint64_t length = rnd_length(gen);
    vec_type test_vec;
    {
      writer_type writer(test_vec);
      for (uint64_t j = 0; j<length; ++j)
        writer << rnd_word(gen);
    }

    std::uniform_int_distribution<> rnd_split(0, length * 64);
    uint64_t split_count = rnd_split(gen);
    std::vector<uint64_t> split_positions;
    split_positions.push_back(0);
    for(uint64_t j = 0; j < split_count; ++j)
      split_positions.push_back(rnd_split(gen));
    split_positions.push_back(length * 64);
    std::sort(split_positions.begin(), split_positions.end());

    std::cout
        << "Running test " << i + 1 << " of " << test_count
        << " (words: " << length
        << ", bits: " << length * 64
        << ", splits: " << split_count << ")" << std::endl;

    vec_type test_vec_restored;
    external_merger merger(test_vec, test_vec_restored);
    for(uint64_t j = 1; j < split_positions.size(); ++j) {
      uint64_t start_pos = split_positions[j - 1];
      uint64_t split_len = split_positions[j] - start_pos;
//      std::cout << "(" << start_pos << "," << split_len << ") ";
      merger.write(split_len, start_pos);
    }
    merger.finish();
//    std::cout << std::endl << std::endl;

    ASSERT_EQ(test_vec_restored.size(), test_vec.size());

    reader_type original_reader(test_vec);
    reader_type restored_reader(test_vec_restored);
//    int k=test_vec.size();
    for(auto original_word : original_reader) {
//      std::cout << k-- << " " << std::flush;
//      ASSERT_EQ(
//        restored_reader.empty(),
//        false
//      );
      if(original_word != *restored_reader) {
        ASSERT_EQ(
            std::bitset<64>(original_word).to_string(),
            std::bitset<64>(*restored_reader).to_string()
        );
      }
      ++restored_reader;

    }
  }


}

/******************************************************************************/
