/*******************************************************************************
 * src/benchmark.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark/algorithm.hpp"
#include "benchmark/file_util.hpp"

// auto filter_wm()

int main(int /*argc*/, char const** /*argv[]*/) {
  auto& algo_list = algorithm_list::get_algorithm_list();
  for (const auto& a : algo_list) {
    a->print_info();
  }
  return 0;
}

/******************************************************************************/
