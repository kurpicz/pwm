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

  uint64_t cur_max_char = (1 << levels);

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

  // The number of 0s at the last level is the number of "even" characters
  if (ContextType::compute_zeros) {
    for (uint64_t i = 0; i < cur_max_char; i += 2) {
      zeros[levels - 1] += last_level_hist[i];
    }
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  for (uint64_t level = levels - 1; level > 0; --level) {
    auto&& this_hist = ctx.hist_at_level(level);
    auto&& next_hist = ctx.hist_at_level(level + 1);
    auto&& borders = ctx.borders_at_shard(level);

    // Update the maximum value of a feasible a bit prefix and update the
    // histogram of the bit prefixes
    cur_max_char >>= 1;
    for (uint64_t i = 0; i < cur_max_char; ++i) {
      this_hist[i] = next_hist[i << 1] + next_hist[(i << 1) + 1];
    }

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    borders[0] = 0;
    for (uint64_t i = 1; i < cur_max_char; ++i) {
      auto const this_rho = ctx.rho(level, i);
      auto const prev_rho = ctx.rho(level, i - 1);

      borders[this_rho] = borders[prev_rho] + this_hist[prev_rho];

      if (ContextType::compute_rho) {
        ctx.set_rho(level - 1, i - 1, prev_rho >> 1);
      }
    }

    // The number of 0s is the position of the first 1 in the previous level
    if (ContextType::compute_zeros) {
      zeros[level - 1] = borders[1];
    }
  }

  reader.rewind();
  // Now we insert the bits with respect to their bit prefixes
  for (uint64_t i = 0; i < size; ++i) {
    const auto cur_char = *reader;
    ++reader;
    for (uint64_t level = levels - 1; level > 0; --level) {
      auto&& borders = ctx.borders_at_shard(level);
      const uint64_t prefix_shift = (levels - level);
      const uint64_t cur_bit_shift = prefix_shift - 1;
      const uint64_t pos = borders[cur_char >> prefix_shift]++;
      bv[level][pos >> 6] |=
          (((cur_char >> cur_bit_shift) & 1ULL) << (63ULL - (pos & 63ULL)));
    }
  }

  ctx.hist_at_level(0)[0] = size;
}

/******************************************************************************/
