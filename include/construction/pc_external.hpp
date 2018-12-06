/*******************************************************************************
 * include/util/pc_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "building_blocks.hpp"

template <typename InputType, typename ContextType>
void pc_in_external(const InputType& text,
                    const uint64_t size,
                    const uint64_t levels,
                    ContextType& ctx) {

  using stxxl_vector_type = InputType;
  using stxxl_reader_type = typename stxxl_vector_type::bufreader_type;

  auto&& zeros = ctx.zeros();
  auto& bv = ctx.bv();

  stxxl_reader_type reader(text);
  auto&& last_level_hist = ctx.hist_at_level(levels);

  // While initializing the histogram, we also compute the first level
  write_bits_wordwise(0, size, bv[0], [&](uint64_t) {
      const auto cur_char = *reader;
      ++reader;
      ++last_level_hist[cur_char];
      return ((cur_char >> (levels - 1)) & 1ULL);
  });

  //~ std::cout << "First level done." << std::endl;

  if constexpr (ContextType::compute_zeros) {
    compute_last_level_zeros(levels, zeros, last_level_hist);
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  bottom_up_compute_hist_borders_optional_zeros_rho(size, levels, ctx);

  reader.rewind();
  // Now we insert the bits with respect to their bit prefixes
  for (uint64_t i = 0; i < size; ++i) {
    const auto cur_char = *reader;
    ++reader;
    for (uint64_t level = levels - 1; level > 0; --level) {
      auto&& borders = ctx.borders_at_level(level);
      const uint64_t prefix_shift = (levels - level);
      const uint64_t cur_bit_shift = prefix_shift - 1;
      const uint64_t pos = borders[cur_char >> prefix_shift]++;
      bv[level][pos >> 6] |=
          (((cur_char >> cur_bit_shift) & 1ULL) << (63ULL - (pos & 63ULL)));
    }
  }
}

/******************************************************************************/
