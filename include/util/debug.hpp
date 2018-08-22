/*******************************************************************************
 * include/util/debug.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <iostream>
#include <string>

#include "util/bit_vectors.hpp"

template<typename WordType, typename bv_t>
std::string bit_string(bv_t const& bv, uint64_t const size) {
  constexpr WordType BITS = (sizeof(WordType) * CHAR_BIT);

  auto s = std::string(size * BITS, '0');
  for (uint64_t bit = 0; bit < size * BITS; bit++) {
    s[bit] += bit_at<WordType>(bv, bit);
  }

  return s;
}

static std::vector<std::vector<uint64_t>> level_sizes(const bit_vectors& bv,
  uint64_t bit_offset, uint64_t bit_length, uint64_t level) {

  if (level == bv.levels()) {
    return {};
  }

  uint64_t zeroes = 0;

  for(uint64_t i = 0; i < bit_length; ++i) {
    if (bit_at(bv[level], bit_offset + i) == 0) {
      ++zeroes;
    }
  }

  uint64_t size_left = zeroes;
  uint64_t size_right = bit_length - zeroes;

  auto sizes_left  = level_sizes(bv, bit_offset, size_left, level + 1);
  auto sizes_right = level_sizes(
    bv, bit_offset + size_left, size_right, level + 1);

  std::vector<std::vector<uint64_t>> r;

  r.push_back({});
  r.back().push_back({bit_length});
  //r.back().push_back({size_right});

  for (uint64_t j = 0; j < sizes_left.size(); j++) {
    r.push_back({});
    auto& v = r.back();

    for(auto& e : sizes_left[j]) {
      v.push_back(e);
    }
    for(auto& e : sizes_right[j]) {
      v.push_back(e);
    }
  }

  return r;
}

static void print_bv(const bit_vectors& bv, uint64_t length) {
  for (uint64_t i = 0; i < bv.size(); i++) {
    std::cout << "   bv["<<i<<"]";

    std::cout << "[";
    for (uint64_t j = 0; j < length; j++) {
      std::cout << uint64_t(bit_at(bv[i], j)) << "";
    }
    std::cout << "]";

    std::cout << "\n";
  }
  std::cout << "\n";
}

static void print_bv_zeros(const bit_vectors& bv,
  const std::vector<uint64_t>& zeros, uint64_t length) {
  for (uint64_t i = 0; i < bv.size(); i++) {
    std::cout << "   bv["<<i<<"]";

    std::cout << "[";
    for (uint64_t j = 0; j < length; j++) {
      std::cout << uint64_t(bit_at(bv[i], j)) << "";
    }
    std::cout << "]";
    std::cout << " zeros[" << i << "] = " << zeros[i];
    std::cout << std::endl;;
  }
  std::cout << std::endl;;
}

[[gnu::unused]] // TODO: C++17 [[maybe_unused]]
static std::string decode_wt(const bit_vectors& bv, uint64_t length) {
  auto ls = level_sizes(bv, 0, length, 0);

  for (auto& v : ls) {
    for (uint64_t i = 1; i < v.size(); i++) {
      v[i] = v[i - 1] + v[i];
    }
    for (uint64_t i = 1; i < v.size(); i++) {
      uint64_t j = v.size() - i;
      v[j] = v[j - 1];
    }
    if (v.size() > 0) {
      v[0] = 0;
    }
  }

  auto r = std::vector<uint8_t>(length);

  for (uint64_t i = 0; i < length; i++) {
    uint8_t value = 0;
    uint64_t j = 0;
    for (uint64_t level = 0; level < bv.levels(); level++) {
      auto& offset = ls[level][j];
      uint8_t bit = bit_at(bv[level], offset);

      value <<= 1;
      value |= bit;

      offset++;
      j = 2 * j + bit;
    }
    r[i] = value;
  }

  return std::string(r.begin(), r.end());
}

[[gnu::unused]] // TODO: C++17 [[maybe_unused]]
static std::string decode_wm(const bit_vectors& bv,
  const std::vector<uint64_t>& zeros, const uint64_t length) {

  if (bv.levels() == 0) {
    return {};
  }

  auto r = std::vector<uint8_t>(length, uint8_t(0));
  auto rtmp = std::vector<uint8_t>(length, uint8_t(0));
  // print_bv_zeros(bv, zeros, length);

  for(uint64_t level = bv.levels() - 1; level > 0; --level) {
    uint64_t offset0 = 0;
    uint64_t offset1 = zeros[level - 1];

    for(uint64_t i = 0; i < length; i++) {
      r[i] |= (bit_at(bv[level], i) << (bv.levels() - level - 1));
    }

    for(uint64_t i = 0; i < length; i++) {
      if(bit_at(bv[level - 1], i) == 0) {
        rtmp[i] = r[offset0];
        offset0++;
      } else {
        rtmp[i] = r[offset1];
        offset1++;
      }
    }

    r.swap(rtmp);
  }

  for(uint64_t i = 0; i < length; i++) {
    r[i] |= bit_at(bv[0], i) << (bv.levels() - 1);
  }

  return std::string(r.begin(), r.end());
}

/******************************************************************************/
