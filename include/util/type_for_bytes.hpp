/*******************************************************************************
 * include/util/type_for_bytes.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "util/debug_assert.hpp"
#include "util/uint_t.hpp"

template <uint8_t BytesPerWord>
struct type_for_bytes {
  type_for_bytes() {
    DCHECK(false); // There must be 1, 2, 4 or 8 bytes per word.
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
struct type_for_bytes<5> {
  using type = uint40_t;
}; // type_for_bytes<4>

template <>
struct type_for_bytes<6> {
  using type = uint48_t;
}; // type_for_bytes<4>

template <>
struct type_for_bytes<8> {
  using type = uint64_t;
}; // type_for_bytes<8>

/******************************************************************************/
