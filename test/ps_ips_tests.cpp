/*******************************************************************************
 * test/construction_test.cpp
 *
 * Copyright (C) 2019 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>

#include <queue>
#include <bitset>
#include <cmath>
#include "util/inplace_partition.hpp"

using vec_t = std::vector<char>;

static void print(const vec_t &vec) {
  for (uint64_t i = 0; i < vec.size(); ++i) {
    std::cout << std::bitset<8>(vec[i]).to_string() << std::endl;
  }
}



TEST(ps_ips, dummy) {

  std::string teststring = "abckjsdfaposadfjlkasd";
  std::vector<char> vec(teststring.begin(), teststring.end());

  auto tesla = ps_ip_sort(vec.data(), vec.size());

  tesla.print();
  tesla.sort(0);
  tesla.print();
  tesla.sort(1);
  tesla.print();
  tesla.sort(2);
  tesla.print();
  tesla.sort(3);
  tesla.print();
  tesla.sort(4);
  tesla.print();
  tesla.sort(5);
  tesla.print();
  tesla.sort(6);
  tesla.print();
  tesla.sort(7);
  tesla.print();

}

/******************************************************************************/
