/*******************************************************************************
 * include/util/ps_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "external_memory/semi_out/inplace_partition.hpp"
#include "external_memory/wavelet_structure_external.hpp"

template <typename AlphabetType, bool is_tree, typename InputType>
void ps_out_external_inplace(
    const InputType& text,
    wavelet_structure_external& result) {

  const auto levels = result.levels();
  const auto size = result.text_size();

  auto& bvs = wavelet_structure_external_writer::bvs(result);
  auto& zeros = wavelet_structure_external_writer::zeros(result);
  auto& hist = wavelet_structure_external_writer::histograms(result);
  std::vector<uint64_t> matrix_sizes(2);

  using result_type = typename std::remove_reference<decltype(bvs)>::type;
  using result_writer_type = typename result_type::bufwriter_type;
  result_writer_type writer(bvs);

  ps_ip_sort partition(text, size, 16ULL * 1024 * 1024);
  partition.start_level();

  // While initializing the histogram, we also compute the first level
  uint64_t cur_pos = 0;
  {
    for (; cur_pos + 64 <= size; cur_pos += 64) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < 64; ++i) {
        word <<= 1;
        const auto character = partition.pop0();
        if constexpr (is_tree) ++hist[levels][character];
        if (((character >> (levels - 1)) & 1ULL)) {
          partition.push1(character);
          word |= 1ULL;
        }
        else
          partition.push0(character);
      }
      writer << word;
    }
    if (size & 63ULL) {
      uint64_t word = 0ULL;
      for (uint64_t i = 0; i < size - cur_pos; ++i) {
        word <<= 1;
        const auto character = partition.pop0();
        if constexpr (is_tree) ++hist[levels][character];
        if (((character >> (levels - 1)) & 1ULL)) {
          partition.push1(character);
          word |= 1ULL;
        }
        else
          partition.push0(character);
      }
      word <<= (64 - (size & 63ULL));
      writer << word;
    }
  }

  partition.finish_level();

  if constexpr (!is_tree) {
    zeros[0] = partition.get_zeros();
  }
  else {
    // compute all histograms
    uint64_t cur_max_char = (1 << levels);
    for (uint64_t level = levels - 1; level > 0; --level) {
      cur_max_char >>= 1;
      for (uint64_t i = 0; i < cur_max_char; ++i) {
        hist[level][i] =
            hist[level + 1][i << 1] + hist[level + 1][(i << 1) + 1];
      }
    }
  }

  std::vector<uint64_t> *sizes = &matrix_sizes;

  // compute all levels (top down)
  for (uint64_t level = 1; level < levels - 1; ++level) {

    if constexpr (!is_tree) {
      matrix_sizes[0] = zeros[level - 1];
      matrix_sizes[1] = size - zeros[level - 1];
    }
    else {
      sizes = &(hist[level]);
    }

    partition.start_level();

    uint64_t i = 0;
    uint64_t word = 0ULL;
    for (uint64_t k = 0; k < sizes->size(); k += 2) {
      const auto node0_size = (*sizes)[k];
      const auto node1_size = (*sizes)[k + 1];

//      std::cout << "Node " << k << " " << node0_size << std::endl;

      for (uint64_t j = 0; j < node0_size; ++j) {
        word <<= 1;
        const auto character = partition.pop0();
        if ((character >> ((levels - 1) - level)) & 1ULL) {
          word |= 1ULL;
          partition.push1(character);
        }
        else {
          partition.push0(character);
        }
        if (PWM_UNLIKELY(++i == 64)) {
          writer << word;
          i = 0;
        }
      }

//      std::cout << "Node " << k + 1 << " " << node1_size << std::endl;

      for (uint64_t j = 0; j < node1_size; ++j) {
        word <<= 1;
        const auto character = partition.pop1();
        if ((character >> ((levels - 1) - level)) & 1ULL) {
          word |= 1ULL;
          partition.push1(character);
        }
        else
          partition.push0(character);
        if (PWM_UNLIKELY(++i == 64)) {
          writer << word;
          i = 0;
        }
      }

    }
    if (i > 0) {
      word <<= (64 - (size & 63ULL));
      writer << word;
    }
    partition.finish_level();

    if constexpr (!is_tree) zeros[level] = partition.get_zeros();
  }

  partition.start_level();
  if constexpr (!is_tree) {
    matrix_sizes[0] = zeros[levels - 2];
    matrix_sizes[1] = size - zeros[levels - 2];
  }
  else {
    sizes = &(hist[levels - 1]);
  }

  uint64_t i = 0;
  uint64_t word = 0ULL;
  for (uint64_t k = 0; k < sizes->size(); k += 2) {
    const auto node0_size = (*sizes)[k];
    const auto node1_size = (*sizes)[k + 1];

    for (uint64_t j = 0; j < node0_size; ++j) {
      word <<= 1;
      word |= (partition.pop0() & 1ULL);
      if (PWM_UNLIKELY(++i == 64)) {
        writer << word;
        i = 0;
      }
    }

    for (uint64_t j = 0; j < node1_size; ++j) {
      word <<= 1;
      word |= (partition.pop1() & 1ULL);
      if (PWM_UNLIKELY(++i == 64)) {
        writer << word;
        i = 0;
      }
    }

  }
  if (i > 0) {
    word <<= (64 - (size & 63ULL));
    writer << word;
  }
}

/******************************************************************************/
