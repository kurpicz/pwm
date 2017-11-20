/*******************************************************************************
 * include/util/histogram.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstdint>
#include <limits>

template <typename AlphabetType>
std::vector<uint64_t> compute_initial_histogram(AlphabetType const* const text,
  const uint64_t size, const uint64_t reduced_sigma = 0) {
  // Compute the histogram
  const uint64_t max_char = std::max(
    static_cast<decltype(reduced_sigma)>(
      std::numeric_limits<AlphabetType>::max() + 1), reduced_sigma);
  std::vector<uint64_t> hist(max_char, 0);
  for (uint64_t pos = 0; pos < size; ++pos) {
    ++hist[text[pos]];
  }
  return hist;
}

/******************************************************************************/
