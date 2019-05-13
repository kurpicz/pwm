/*******************************************************************************
 * include/util/pc_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <omp.h>
#include <cmath>
#include "construction/building_blocks.hpp"
#include "util/cacheline.hpp"

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


template <typename InputType, typename ContextType,
          uint64_t block_bytes = 32ULL * 1024 * 1024,
          typename stats_type>
void pc_in_external_parallel(const InputType& text,
                    const uint64_t size,
                    const uint64_t levels,
                    ContextType& ctx,
                    stats_type &stats) {

  using value_type = typename InputType::value_type;
  using stxxl_vector_type = InputType;
  using stxxl_reader_type = typename stxxl_vector_type::bufreader_type;

  auto&& zeros = ctx.zeros();
  auto& bv = ctx.bv();

  omp_lock_t mutex_read_text;
  omp_lock_t mutex_write_hist;
  omp_init_lock(&mutex_read_text);
  omp_init_lock(&mutex_write_hist);
  uint64_t threads;
  #pragma omp parallel
  threads = (uint64_t)omp_get_num_threads();


  const uint64_t total_chars =
      std::max(threads << 1, block_bytes * threads / sizeof(value_type));
  const uint64_t total_bytes = total_chars * sizeof(value_type);

  value_type * raw_buffer = (value_type *)malloc(total_bytes);

  const uint64_t small_size = total_chars / threads;
  const uint64_t small_blocks = (size + small_size - 1) / small_size;
  const uint64_t last_small_size = size - (small_blocks - 1) * small_size;

  const uint64_t big_size = total_chars >> 1;
  const uint64_t big_blocks = (size + big_size - 1) / big_size;
  const uint64_t last_big_size = size - (big_blocks - 1) * big_size;

  uint64_t next_block = 0;
  stxxl_reader_type reader(text);
  auto&& ll_hist = ctx.hist_at_level(levels);
  const uint64_t ll_hist_size = 1ULL << levels;


//  std::cout << "Starting histogram phase.\n"
//            << "T:  " << total_chars << " / " << total_bytes << "\n"
//            << "B:  " << big_blocks << "\n"
//            << "S:  " << big_size << "\n"
//            << "LS: " << last_big_size << "\n"
//            << "b:  " << small_blocks << "\n"
//            << "s:  " << small_size << "\n"
//            << "ls: " << last_small_size << std::endl;

  stats.phase("scan1");

  #pragma omp parallel
  {
    const uint64_t omp_rank = omp_get_thread_num();

    value_type * small_buffer = &(raw_buffer[omp_rank * small_size]);

    while (true) {
      omp_set_lock(&mutex_read_text);
      if (next_block == small_blocks) {
//        std::cout << "Thread " << omp_rank << " done." << std::endl;
        omp_unset_lock(&mutex_read_text);
        break;
      }
      ++next_block;
      const uint64_t s =
          (next_block == small_blocks) ? last_small_size : small_size;

//      std::cout << "Thread " << omp_rank
//                << ", read block " << next_block - 1
//                << " (" << s << ")" << std::endl;

      for (uint64_t i = 0; i < s; ++i) {
        reader >> small_buffer[i];
      }

//      std::cout << "Thread " << omp_rank
//                << ", done read block " << next_block << std::endl;
      omp_unset_lock(&mutex_read_text);

      std::vector<uint64_t> local_ll_hist(ll_hist_size, 0);
      for (uint64_t i = 0; i < s; ++i) {
        ++local_ll_hist[small_buffer[i]];
      }

      omp_set_lock(&mutex_write_hist);
      for (uint64_t i = 0; i < ll_hist_size; ++i) {
        ll_hist[i] += + local_ll_hist[i];
      }
      omp_unset_lock(&mutex_write_hist);
    }
  }

  if constexpr (ContextType::compute_zeros) {
    compute_last_level_zeros(levels, zeros, ll_hist);
  }

  // Now we compute the WM bottom-up, i.e., the last level first
  bottom_up_compute_hist_borders_optional_zeros_rho(size, levels, ctx);

//  std::cout << "Histogram phase complete." << std::endl;

  value_type * big_buffer1 = raw_buffer;
  value_type * big_buffer2 = &(raw_buffer[big_size]);

  // make sure that the borders of different levels
  // are in different cache lines at all times
  const uint64_t cache_line_bytes = get_cache_line_size();
  const uint64_t cache_line_words = (cache_line_bytes + 7) / 8;
  std::vector<std::vector<uint64_t>> all_borders(levels);
  for (uint64_t l = 0; l < levels; ++l) {
    const uint64_t size = 1ULL << l;
    all_borders[l].resize(size + cache_line_words, 0ULL);
    auto &&borders = ctx.borders_at_level(l);
    for (uint64_t i = 0; i < size; ++i) {
      all_borders[l][i] = borders[i];
    }
  }

  stats.phase("scan2");

  reader.rewind();
  uint64_t s = (big_blocks == 1) ? last_big_size : big_size;
  for (uint64_t i = 0; i < s; ++i) {
    reader >> big_buffer2[i];
  }

  for (uint64_t b = 0; b < big_blocks; ++b) {
    std::swap(big_buffer1, big_buffer2);
    const uint64_t next_s =
        ((b + 2 < big_blocks) ? big_size :
        ((b + 1 < big_blocks) ? last_big_size : 0));

    #pragma omp parallel for schedule(nonmonotonic : dynamic, 1)
    for (uint64_t l = 0; l <= levels; ++l) {
      if (PWM_UNLIKELY(l == 0)) {
        for (uint64_t i = 0; i < next_s; ++i) {
          reader >> big_buffer2[i];
        }
      }
      else {
        const uint64_t level = l - 1;
        auto& borders = all_borders[level];
        const uint64_t prefix_shift = (levels - level);
        const uint64_t cur_bit_shift = prefix_shift - 1;
        auto lbv = bv[level];

        for (uint64_t i = 0; i < s; ++i) {
          const auto cur_char = big_buffer1[i];
          const uint64_t pos = borders[cur_char >> prefix_shift]++;
          lbv[pos >> 6] |=
              (((cur_char >> cur_bit_shift) & 1ULL) << (63ULL - (pos & 63ULL)));
        }
      }
    }
    s = next_s;
  }
  delete raw_buffer;
}

/******************************************************************************/
