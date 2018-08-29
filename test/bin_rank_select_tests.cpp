/*******************************************************************************
 * test/bin_rank_select.cpp
 *
 * Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>

#include <random>

#include "arrays/bit_vectors.hpp"
#include "queries/bin_rank_select_popcnt.hpp"
#include "util/common.hpp" 

TEST(bin_rank_select_popcnt, rank) {
  const size_t bit_vector_size = 1024*1024*128;

  bit_vectors bv(1, bit_vector_size);

  std::mt19937 gen(44227);
  std::uniform_int_distribution<uint64_t> dis(0);
  for (size_t i = 0; i < word_size(bit_vector_size); ++i) {
    *(bv[0] + i) = dis(gen);
  }

  bin_rank_select_popcnt rank_support(bv[0], word_size(bit_vector_size));

  size_t one_count = 0;
  for (size_t i = 0; i < bit_vector_size; ++i) {
    ASSERT_EQ(rank_support.rank0(i), i - one_count);
    ASSERT_EQ(rank_support.rank1(i), one_count);
    if (bit_at(bv[0], i)) { one_count++; }
  }
}

/******************************************************************************/
