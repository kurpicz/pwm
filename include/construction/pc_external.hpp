/*******************************************************************************
 * include/util/pc_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/
 
#pragma once

template <typename InputType, typename ContextType>
void pc_in_external(const InputType& text, const uint64_t size, const uint64_t levels,
  ContextType& ctx) {

  using stxxl_vector_type = InputType;
  using stxxl_reader_type = typename stxxl_vector_type::bufreader_type;
  using borders_type = typename std::remove_reference<decltype(ctx.borders())>::type;
  uint64_t borders_size = ctx.borders().size();

  uint64_t cur_max_char = (1 << levels);
  uint64_t cur_pos = 0;
  
  auto& zeros = ctx.zeros();
  auto& bv = ctx.bv();
  std::vector<borders_type> borders_v(levels);
  borders_type * borders = borders_v.data();
    
  stxxl_reader_type reader(text);
  // While initializing the histogram, we also compute the first level
  for (; cur_pos + 64 <= size; cur_pos += 64) {
    uint64_t word = 0ULL;
    for (uint64_t i = 0; i < 64; ++i) {
      const auto cur_char = *reader;
      ++reader;
      ++ctx.hist(levels, cur_char);
      word <<= 1;
      word |= ((cur_char >> (levels - 1)) & 1ULL);
    }
    bv[0][cur_pos >> 6] = word;
  }
  if (size & 63ULL) {
    uint64_t word = 0ULL;
    for (uint64_t i = 0; i < size - cur_pos; ++i) {
      const auto cur_char = *reader;
      ++reader;
      ++ctx.hist(levels, cur_char);
      word <<= 1;
      word |= ((cur_char >> (levels - 1)) & 1ULL);
    }
    word <<= (64 - (size & 63ULL));
    bv[0][size >> 6] = word;
  }

  //~ std::cout << "First level done." << std::endl;

  // The number of 0s at the last level is the number of "even" characters
  if (ContextType::compute_zeros) {
    for (uint64_t i = 0; i < cur_max_char; i += 2) {
      zeros[levels - 1] += ctx.hist(levels, i);
    }
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  for (uint64_t level = levels - 1; level > 0; --level) {

    // Update the maximum value of a feasible a bit prefix and update the
    // histogram of the bit prefixes
    cur_max_char >>= 1;
    for (uint64_t i = 0; i < cur_max_char; ++i) {
      ctx.hist(level, i)
        = ctx.hist(level + 1, i << 1)
        + ctx.hist(level + 1, (i << 1) + 1);
    }

    borders[level].resize(borders_size);
    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    borders[level][0] = 0;
    for (uint64_t i = 1; i < cur_max_char; ++i) {
      auto const prev_rho = ctx.rho(level, i - 1);

      borders[level][ctx.rho(level, i)] =
        borders[level][prev_rho] + ctx.hist(level, prev_rho);

      if (ContextType::compute_rho)  {
        ctx.set_rho(level - 1, i - 1, prev_rho >> 1);
      }
    }

    // The number of 0s is the position of the first 1 in the previous level
    if (ContextType::compute_zeros) {
      zeros[level - 1] = borders[level][1];
    }
  }

  reader.rewind();
  // Now we insert the bits with respect to their bit prefixes
  for (uint64_t i = 0; i < size; ++i) {
    const auto cur_char = *reader;
    ++reader;
    for (uint64_t level = levels - 1; level > 0; --level) {
      const uint64_t prefix_shift = (levels - level);
      const uint64_t cur_bit_shift = prefix_shift - 1;
      const uint64_t pos = borders[level][cur_char >> prefix_shift]++;
      bv[level][pos >> 6] |= (((cur_char >> cur_bit_shift) & 1ULL)
        << (63ULL - (pos & 63ULL)));
    }
  }

  if (levels > 1) { // TODO check condition
    ctx.hist(0, 0) = ctx.hist(1, 0) + ctx.hist(1, 1);
  }
}

/******************************************************************************/
