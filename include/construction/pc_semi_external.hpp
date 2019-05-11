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

template <typename InputType, typename ContextType>
void pc_in_external_parallel(InputType& text_raw,
                             const uint64_t size,
                             const uint64_t levels,
                             ContextType& ctx) {

  using ctx_t = ContextType;
  using char_t = typename InputType::value_type;
  using stxxl_vector_type = InputType;
  using stxxl_reader_type = typename stxxl_vector_type::bufreader_type;
  using stxxl_const_iter_type = typename stxxl_vector_type::const_iterator;

  constexpr static uint64_t block_bytes = 2ULL * 1024 *1024;
  constexpr static uint64_t block_chars = block_bytes / sizeof(char_t);
  constexpr static uint64_t page_blocks = 4ULL;
  constexpr static uint64_t page_bytes = block_bytes * page_blocks;
  constexpr static uint64_t page_chars = block_chars * page_blocks;

  auto& bv = ctx.bv();
  auto& zeros = ctx.zeros();

  struct access_type {
    const InputType& text_;
    const uint64_t lbound_;
    const uint64_t rbound_;

    stxxl_const_iter_type it_;
    char_t * page_data_;
    uint64_t page_start_;
    uint64_t page_current_;

    access_type(const InputType& text,
                const uint64_t lbound,
                const uint64_t rbound)
        : text_(text),
          lbound_(lbound),
          rbound_(rbound),
          it_(text_.cbegin() + lbound),
          page_start_(0),
          page_current_(0) {
      page_data_ = static_cast<char_t *>(malloc(page_bytes));
      load_page();
    }

    access_type(const InputType& text) : access_type(text, 0, text.size()) {}

    void operator >> (char_t& v) {
      v = page_data_[page_current_++];
      if(PWM_UNLIKELY(page_current_ == page_chars)) {
        page_start_ += page_chars;
        page_current_ = 0;
        load_page();
      }
    }

    void load_page() {
      const uint64_t count = std::min(rbound_ - page_start_, page_chars);
      #pragma omp critical
      {
        for (uint64_t i = 0; i < count; ++i) {
          page_data_[i] = *it_;
          ++it_;
        }
      }
    }

    ~access_type() {
      delete page_data_;
    }
  };

  std::vector<uint64_t> initial_hist((1ULL << levels) * levels, 0);

  #pragma omp parallel num_threads(levels)
  {
    const uint64_t omp_rank = uint64_t(omp_get_thread_num());
    const uint64_t omp_size = uint64_t(omp_get_num_threads());
    const uint64_t alphabet_size = (1 << levels);
    const bool is_tail = (omp_rank + 1 == omp_size);

    auto *const initial_hist_ptr =
        initial_hist.data() + (alphabet_size * omp_rank);


    const uint64_t slice_size = ((size / omp_size) >> 6) << 6;
    const uint64_t lbound = omp_rank * slice_size;
    const uint64_t rbound = is_tail ? (size - 63) : (lbound + slice_size);

    {
      access_type slice_reader(text_raw, lbound, rbound);
      char_t character;
      for (uint64_t cur_pos = lbound; cur_pos < rbound; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < 64; ++i) {
          slice_reader >> character;
          ++initial_hist_ptr[character];
          word <<= 1;
          word |= ((character >> (levels - 1)) & 1ULL);
        }
        bv[0][cur_pos >> 6] = word;
      }

      if ((size & 63ULL) && is_tail) {
        uint64_t word = 0ULL;
        for (uint64_t i = 0; i < (size & 63ULL); ++i) {
          slice_reader >> character;
          ++initial_hist_ptr[character];
          word <<= 1;
          word |= ((character >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        bv[0][size >> 6] = word;
      }
    }

    #pragma omp barrier

    // Compute the historam with respect to the local slices of the text
    #pragma omp for
    for (uint64_t i = 0; i < alphabet_size; ++i) {
      for (uint64_t rank = 0; rank < omp_size; ++rank) {
        ctx.hist(levels, i) +=
            *(initial_hist.data() + (alphabet_size * rank) + i);
      }
    }

    #pragma omp single
    {
      if constexpr (ctx_t::compute_zeros) {
        // The number of 0s at the last level is the number of "even"
        // characters
        for (uint64_t i = 0; i < alphabet_size; i += 2) {
          zeros[levels - 1] += ctx.hist(levels, i);
        }
      }
    }

    // Compute the histogram for each level of the wavelet structure
    #pragma omp for
    for (uint64_t level = 1; level < levels; ++level) {
      const uint64_t local_alphabet_size = (1 << level);
      const uint64_t requierd_characters = (1 << (levels - level));
      for (uint64_t i = 0; i < local_alphabet_size; ++i) {
        for (uint64_t j = 0; j < requierd_characters; ++j) {
          ctx.hist(level, i) +=
              ctx.hist(levels, (i * requierd_characters) + j);
        }
      }
    }
  }


  char_t * page_load = static_cast<char_t *>(malloc(page_bytes));
  char_t * page_use = static_cast<char_t *>(malloc(page_bytes));
  uint64_t page_offset = 0;
  uint64_t page_size = 0;
  const uint64_t page_count = (text_raw.size() + page_chars - 1) / page_chars;

  stxxl_reader_type reader(text_raw);

  std::vector<uint64_t> prefix_shift(levels);
  std::vector<uint64_t> cur_bit_shift(levels);
  std::vector<std::vector<uint64_t>> borders(levels);

    // Now we compute the wavelet structure bottom-up, i.e., the last level
    // first
  #pragma omp parallel num_threads(levels)
  {
    const uint64_t level = uint64_t(omp_get_thread_num());
    const uint64_t local_alphabet_size = (1 << level);
    prefix_shift[level] = (levels - level);
    cur_bit_shift[level] = prefix_shift[level] - 1;
    borders[level].resize(local_alphabet_size);

    if (level == 0) {
      page_size = std::min(uint64_t(text_raw.size()) - page_offset, page_chars);
      page_offset += page_size;
      for (uint64_t i = 0; i < page_size; ++i) reader >> page_use[i];
    }
    else {
      borders[level][0] = 0;
      for (uint64_t i = 1; i < local_alphabet_size; ++i) {
        const auto prev_rho = ctx.rho(level, i - 1);
        borders[level][ctx.rho(level, i)] =
            borders[level][prev_rho] + ctx.hist(level, prev_rho);
      }
      // The number of 0s is the position of the first 1 in the previous level
      if constexpr (ctx_t::compute_zeros) {
        zeros[level - 1] = borders[level][1];
      }
    }
  }

  for (uint64_t page = 0; page < page_count; ++page) {
    const uint64_t next_page_size = std::min(uint64_t(text_raw.size()) - page_offset, page_chars);

    #pragma omp parallel num_threads(levels)
    {
      const uint64_t level = uint64_t(omp_get_thread_num());

      if (level == 0) {
        for (uint64_t i = 0; i < next_page_size; ++i) reader >> page_load[i];
      }
      else {
        const uint64_t prefix_shiftl = prefix_shift[level];
        const uint64_t cur_bit_shiftl = cur_bit_shift[level];
        std::vector<uint64_t> &bordersl = borders[level];
        for (uint64_t i = 0; i < page_size; ++i) {
          const uint64_t pos = bordersl[page_use[i] >> prefix_shiftl]++;
          bv[level][pos >> 6] |=
              (((page_use[i] >> cur_bit_shiftl) & 1ULL) << (63ULL - (pos & 63ULL)));
        }
      }
    }

    page_size = next_page_size;
    page_offset += next_page_size;
    std::swap(page_load, page_use);
  }

  delete page_load;
  delete page_use;
}

/******************************************************************************/
