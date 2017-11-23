/*******************************************************************************
 * include/util/file_util.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <fstream>
#include <vector>

#include "type_for_bytes.hpp"

template <uint8_t BytesPerWord>
static std::vector<typename type_for_bytes<BytesPerWord>::type> file_to_vector(
  const std::string& file_name) {
  std::ifstream stream(file_name.c_str(), std::ios::in | std::ios::binary);
  stream.seekg(0, std::ios::end);
  uint64_t size = stream.tellg() / BytesPerWord;
  stream.seekg(0);
  std::vector<typename type_for_bytes<BytesPerWord>::type> result(size);
  stream.read(reinterpret_cast<char*>(result.data()), size);
  stream.close();
  return result;
}

template <typename Type>
static void vector_to_file(const std::vector<Type>& vec,
  const std::string& file_name) {

  std::ofstream stream(file_name.c_str(), std::ios::binary | std::ios::out);
  stream.write(reinterpret_cast<char*>(vec.data()), sizeof(Type) * vec.size());
  stream.close();
}

/******************************************************************************/
