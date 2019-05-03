/*******************************************************************************
 * include/wx_ps_fe.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include <thread>
#include <sstream>
#include <omp.h>

#include "arrays/memory_types.hpp"
#include "construction/pc_dd_fe/merge_external.hpp"
#include "construction/pc_dd_fe/pc_partial.hpp"
#include "construction/pc_dd_fe/ctx_partial.hpp"
#include "construction/wavelet_structure.hpp"
#include "construction/wavelet_structure_external.hpp"
#include "util/permutation.hpp"

#include "wx_base.hpp"
#include "wx_dd_pc.hpp"
#include "wx_pc.hpp"

#define DDE2_VERBOSE if constexpr (true) atomic_out()


template <typename AlphabetType, bool is_tree_, uint64_t bytes_memory_ = 1024 * 1024 * 1024>
class wx_dd_pc_fe2 : public wx_in_out_external<true, true, true> {

  class atomic_out {
    std::ostringstream stream;
  public:
    template <typename in_type>
    atomic_out& operator<<(const in_type& in) {
      stream << in;
      return *this;
    }
    ~atomic_out() {
      std::cout << stream.str() << std::flush;
    }
  };

  static uint64_t get_omp_size() {
    int result;
    #pragma omp parallel
    result = omp_get_num_threads();
    return (uint64_t)result;
  }

  template <typename InputType>
  struct dd_ctx {
    using char_type = typename InputType::value_type;
    using result_type = typename wavelet_structure_external::bv_type;
    using result_writer_type = typename result_type::bufwriter_type;
    using text_reader_type = typename InputType::bufreader_type;

    const uint64_t omp_size_;

    const InputType& text_;
    const uint64_t size_;
    const uint64_t levels_;

    const uint64_t input_bytes_per_block_;
    const uint64_t input_chars_per_block_;
    const uint64_t number_of_blocks_;
    const uint64_t input_chars_last_block_;

    std::vector<pow2_array> block_hists_;
    std::deque<char_type *> text_blocks_;

    void run() {
      text_reader_type reader(text_);
      block_hists_.reserve(number_of_blocks_);

      char_type** read = new char_type*[number_of_blocks_ + omp_size_];
      ctx_partial** compute = new ctx_partial*[number_of_blocks_ + omp_size_];
      bool* write = new bool[number_of_blocks_ + omp_size_];

      #pragma omp parallel
      {
      #pragma omp single
      {
      for (uint64_t b = 0; b < number_of_blocks_; ++b) {

        #pragma omp task depend(in: read[b + omp_size_ - 1], compute[b]) depend(out: read[b + omp_size_])
        {
          const uint64_t block_size = ((b + 1) == number_of_blocks_) ?
                                      input_chars_last_block_ : input_chars_per_block_;

          char_type * text_block = text_blocks_.front();
          text_blocks_.pop_front();
          text_blocks_.push_back(text_block);
          for (uint64_t i = 0; i < block_size; ++i) {
            reader >> text_block[i];
          }
          read[b + omp_size_] = text_block;

          DDE2_VERBOSE << "Read " << b << " (bs " << block_size << ")\n";
        }

        #pragma omp task depend(in: read[b + omp_size_], write[b]), depend(out: compute[b + omp_size_])
        {
          const uint64_t block_size = ((b + 1) == number_of_blocks_) ?
                                      input_chars_last_block_ : input_chars_per_block_;

          char_type * text_block = read[b + omp_size_];
          ctx_partial * ctx_block = new ctx_partial(block_size, levels_);
          pc_partial(text_block, ctx_block);
          compute[b + omp_size_] = ctx_block;

          DDE2_VERBOSE << "Compute " << b << " " << text_block[1] << " (bs " << block_size << ")\n";
        }

        #pragma omp task depend(in: compute[b + omp_size_], write[b + omp_size_ - 1]) depend(out: write[b + omp_size_])
        {
          const uint64_t block_size = ((b + 1) == number_of_blocks_) ?
                                      input_chars_last_block_ : input_chars_per_block_;
          ctx_partial * ctx_block = compute[b + omp_size_];
          block_hists_.push_back(ctx_block->extract_final_hist());

          delete compute[b + omp_size_];
          DDE2_VERBOSE << "Write " << b << " (bs " << block_size << ")\n";
        }
      }
      }
      }

      delete read;
      delete compute;
      delete write;

    }

    dd_ctx(const InputType &text, const uint64_t size, const uint64_t levels)
        : omp_size_(get_omp_size()),
          text_(text), size_(size), levels_(levels),
          input_bytes_per_block_((bytes_memory_ / omp_size_) >> 1),
          input_chars_per_block_(input_bytes_per_block_ / sizeof(char_type)),
          number_of_blocks_((size + input_chars_per_block_ - 1) / input_chars_per_block_),
          input_chars_last_block_(size_ - input_chars_per_block_ * (number_of_blocks_ - 1)) {

      for (uint64_t i = 0; i < omp_size_; ++i) {
        text_blocks_.push_back((char_type*)malloc(input_chars_per_block_ * sizeof(char_type)));
      }

      DDE2_VERBOSE
        << "Created context for wx_dd_fe "
        << "[blocks: " << number_of_blocks_ << ", "
        << "block_chars: " << input_chars_per_block_ << ", "
        << "last_block_chars: " << input_chars_last_block_  << ", "
        << "omp_size: " << omp_size_
        << "].\n";
    }

    ~dd_ctx() {
      while(!text_blocks_.empty()) {
        delete text_blocks_.front();
        text_blocks_.pop_front();
      }
    }
  };

public:
  static constexpr bool is_parallel = true;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  template <typename InputType, typename stats_type>
  static wavelet_structure_external
  compute(const InputType& text,
          const uint64_t size,
          const uint64_t levels,
          stats_type& stats) {
    static_assert(sizeof(typename InputType::value_type) == sizeof(AlphabetType));

    stats.phase("dd");

    // create empty result
    std::ostringstream name;
    name << "w" << (is_tree_ ? "t" : "m") << "_dd_fe";
    auto result =
        wavelet_structure_external_factory(is_tree_).
            histograms(is_tree_).zeros(!is_tree_).
            construct(size, levels, name.str(), 0);
    if(size == 0) return result;

    wavelet_structure_external_writer::bvs(result).clear();

    using ctx_t = dd_ctx<InputType>;
    ctx_t ctx(text, size, levels);

    ctx.run();

    return result;
  }
}; // class wx_ps

/******************************************************************************/
