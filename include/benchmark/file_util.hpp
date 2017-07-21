/*******************************************************************************
 * include/file_util.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef FILE_UTIL_HEADER
#define FILE_UTIL_HEADER

#include <cassert>
#include <fstream>
#include <vector>

template <uint8_t BytesPerWord>
struct type_for_bytes {
  type_for_bytes() {
    assert(false); // There must be 1, 2, 4 or 8 bytes per word.
  }
}; // type_for_bytes<uint8_t>

template <>
struct type_for_bytes<1> {
  using type = uint8_t;
}; // type_for_bytes<1>

template <>
struct type_for_bytes<2> {
  using type = uint16_t;
}; // type_for_bytes<2>

template <>
struct type_for_bytes<4> {
  using type = uint32_t;
}; // type_for_bytes<4>

template <>
struct type_for_bytes<8> {
  using type = uint64_t;
}; // type_for_bytes<8>

template <uint8_t BytesPerWord>
std::vector<typename type_for_bytes<BytesPerWord>::type> file_to_vector(
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

#endif // FILE_UTIL_HEADER

/******************************************************************************/
