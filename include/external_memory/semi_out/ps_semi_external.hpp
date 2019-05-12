/*******************************************************************************
 * include/util/ps_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "external_memory/wavelet_structure_external.hpp"

template <typename AlphabetType, bool is_tree, typename InputType>
void ps_out_external(
    const InputType& text,
    wavelet_structure_external& result) {
  auto levels = result.levels();
  auto size = result.text_size();

  auto& bvs = wavelet_structure_external_writer::bvs(result);
  auto& zeros = wavelet_structure_external_writer::zeros(result);
  auto& hist = wavelet_structure_external_writer::histograms(result);
  auto const& level_offsets = result.level_offsets();

  std::vector<uint64_t> borders(1ULL << levels, 0);
  auto rho = rho_dispatch<is_tree>::create(levels);

  auto sorted_text_vec = std::vector<AlphabetType>(size);
  AlphabetType* sorted_text = sorted_text_vec.data();

  uint64_t cur_max_char = (1 << levels);

  using result_type = typename std::remove_reference<decltype(bvs)>::type;
  using result_writer_type = typename result_type::bufwriter_type;

  // While initializing the histogram, we also compute the first level
  uint64_t cur_pos = 0;
  {
    result_writer_type writer(bvs);
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < 64; ++i) {
        ++hist[levels][text[cur_pos + i]];
        word <<= 1;
        word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      }
      writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < size - cur_pos; ++i) {
        ++hist[levels][text[cur_pos + i]];
        word <<= 1;
        word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
      }
      word <<= (64 - (size & 63ULL));
      writer << word;
    }
  }

  // The number of 0s at the last level is the number of "even" characters
  if constexpr (!is_tree) {
    for (uint64_t i = 0; i < cur_max_char; i += 2) {
      zeros[levels - 1] += hist[levels][i];
    }
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  for (uint64_t level = levels - 1; level > 0; --level) {
    // Update the maximum value of a feasible a bit prefix and update the
    // histogram of the bit prefixes
    cur_max_char >>= 1;
    for (uint64_t i = 0; i < cur_max_char; ++i) {
      hist[level][i] =
          hist[level + 1][i << 1] + hist[level + 1][(i << 1) + 1];
    }

    // Compute the starting positions of characters with respect to their
    // bit prefixes and the bit-reversal permutation
    borders[0] = 0;
    for (uint64_t i = 1; i < cur_max_char; ++i) {
      auto const prev_rho = rho(level, i - 1);

      borders[rho(level, i)] =
          borders[prev_rho] + hist[level][prev_rho];
    }

    // The number of 0s is the position of the first 1 in the previous level
    if constexpr (!is_tree) {
      zeros[level - 1] = borders[1];
    }

    // Now we sort the text utilizing counting sort and the starting positions
    // that we have computed before
    for (uint64_t i = 0; i < size; ++i) {
      const auto cur_char = text[i];
      sorted_text[borders[cur_char >> (levels - level)]++] = cur_char;
    }

    {
      result_writer_type writer(bvs.begin() + level_offsets[level]);
      // Since we have sorted the text, we can simply scan it from left to right
      // and for the character at position $i$ we set the $i$-th bit in the
      // bit vector accordingly
      for (cur_pos = 0; cur_pos + 63 < size; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < 64; ++i) {
          word <<= 1;
          word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
        }
        writer << word;
      }
      if (size & 63ULL) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < size - cur_pos; ++i) {
          word <<= 1;
          word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        writer << word;
      }
    }
  }
}

/******************************************************************************/
