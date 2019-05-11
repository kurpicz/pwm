/*******************************************************************************
 * include/util/meminfo.hpp
 *
 * Copyright (C) 2019 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <fstream>

uint64_t get_meminfo(const std::string row) {
  const std::string header = row + ":";
  std::string token;
  std::ifstream file("/proc/meminfo");
  while(file >> token) {
    if(token == header) {
      uint64_t res;
      if(file >> res) {
        return res;
      } else {
        return 0;
      }
    }
  }
  return 0; // nothing found
}


uint64_t get_mem_total() {
  return get_meminfo("MemTotal");
}

uint64_t get_mem_available() {
  return get_meminfo("MemAvailable");
}

/******************************************************************************/
