/*******************************************************************************
 * include/util/file_util.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <fstream>
#include <vector>

#include "util/type_for_bytes.hpp"

template <uint8_t BytesPerWord>
static std::vector<typename type_for_bytes<BytesPerWord>::type>
file_to_vector(const std::string& file_name, const size_t prefix_size = 0) {
  std::ifstream stream(file_name.c_str(), std::ios::in | std::ios::binary);

  if (!stream) {
    std::cerr << "File " << file_name << " not found\n";
    exit(1);
  }

  stream.seekg(0, std::ios::end);
  uint64_t size = stream.tellg() / BytesPerWord;
  if (prefix_size > 0) {
    size = std::min(prefix_size, size);
  }
  stream.seekg(0);
  std::vector<typename type_for_bytes<BytesPerWord>::type> result(size);
  stream.read(reinterpret_cast<char*>(result.data()), size);
  stream.close();
  return result;
}

/******************************************************************************/
