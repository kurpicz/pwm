/*******************************************************************************
 * include/construction/pc_dd_fe/wx_dd_pc_fe2.hpp
 *
 * Copyright (C) 2019 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <bitset>
#include <vector>
#include <thread>
#include <sstream>
#include <omp.h>

#include "arrays/memory_types.hpp"
#include "construction/pc_dd_fe/external_merge.hpp"
#include "construction/pc_dd_fe/pc_partial.hpp"
#include "construction/pc_dd_fe/ctx_partial.hpp"
#include "construction/wavelet_structure.hpp"
#include "construction/wavelet_structure_external.hpp"
#include "util/permutation.hpp"
#include "util/print.hpp"

#include "wx_base.hpp"
#include "wx_dd_pc_fe.hpp"
#include "wx_pc.hpp"

#define DDE2_VERBOSE if constexpr (false) atomic_out()


template <typename AlphabetType, bool is_tree_, uint64_t bytes_memory_ = 1024 * 1024 * 1024, bool rw_simultaneously = false>
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
    using result_reader_type = typename result_type::bufreader_type;
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

    std::vector<std::vector<uint64_t>> block_hists_;

    result_type unmerged_result;

    uint64_t get_chars_per_block(const uint64_t block) {
      return ((block + 1) == number_of_blocks_) ?
             input_chars_last_block_ : input_chars_per_block_;
    }

    void run() {
      text_reader_type reader(text_);
      result_writer_type writer(unmerged_result);

      std::deque<char_type *> text_blocks;
      for (uint64_t i = 0; i < omp_size_; ++i) {
        text_blocks.push_back((char_type*)malloc(input_chars_per_block_ * sizeof(char_type)));
      }

      char_type** read = new char_type*[number_of_blocks_ + omp_size_];
      ctx_partial** compute = new ctx_partial*[number_of_blocks_ + omp_size_];
      bool* write = new bool[number_of_blocks_ + omp_size_];

      omp_lock_t rw_lock;
      if constexpr (!rw_simultaneously)
        omp_init_lock(&rw_lock);

      #pragma omp parallel
      {
      #pragma omp single
      {
      for (uint64_t b = 0; b < number_of_blocks_; ++b) {

        #pragma omp task depend(in: read[b + omp_size_ - 1], compute[b]) depend(out: read[b + omp_size_])
        {
          if constexpr (!rw_simultaneously) omp_set_lock(&rw_lock);
          const uint64_t block_size = get_chars_per_block(b);
          char_type * text_block = text_blocks.front();
          text_blocks.pop_front();
          text_blocks.push_back(text_block);
          for (uint64_t i = 0; i < block_size; ++i) {
            reader >> text_block[i];
          }
          read[b + omp_size_] = text_block;

          DDE2_VERBOSE << "Read " << b << " (bs " << block_size << ")\n";
          if constexpr (!rw_simultaneously) omp_unset_lock(&rw_lock);
        }

        #pragma omp task depend(in: read[b + omp_size_], write[b]), depend(out: compute[b + omp_size_])
        {
          const uint64_t block_size = get_chars_per_block(b);
          char_type * text_block = read[b + omp_size_];
          DDE2_VERBOSE << "START Compute " << b << " " << text_block[1] << " (bs " << block_size << ")\n";
          ctx_partial * ctx_block = new ctx_partial(block_size, levels_);
          pc_partial(text_block, ctx_block);
          compute[b + omp_size_] = ctx_block;
          DDE2_VERBOSE << "DONE Compute " << b << " " << text_block[1] << " (bs " << block_size << ")\n";
        }

        #pragma omp task depend(in: compute[b + omp_size_], write[b + omp_size_ - 1]) depend(out: write[b + omp_size_])
        {
          if constexpr (!rw_simultaneously) omp_set_lock(&rw_lock);
          DDE2_VERBOSE << "START Write " << b << " (bs " << get_chars_per_block(b) << ")\n";
          ctx_partial * ctx_block = compute[b + omp_size_];
          block_hists_.emplace_back(ctx_block->flat_hist());
          const auto &block_bvs = ctx_block->bv();
          const auto &level_data_sizes = ctx_block->level_data_sizes();
          DDE2_VERBOSE << "Levels " << levels_ << ", ds " << ctx_block->data_size() << "\n";
          for (uint64_t l = 0; l < levels_; ++l) {
            const auto level_data_size = level_data_sizes[l];
            DDE2_VERBOSE << "    Level " << l << ", ds " << level_data_size << "\n";
            for (uint64_t w = 0; w < level_data_size; ++w) {
              writer << block_bvs[l][w];
              DDE2_VERBOSE << std::bitset<64>(block_bvs[l][w]).to_string() << "\n";
            }
          }
          delete ctx_block;
          DDE2_VERBOSE << "DONE Write " << b << " (bs " << get_chars_per_block(b) << ")\n";
          if constexpr (!rw_simultaneously) omp_unset_lock(&rw_lock);
        }
      }
      }
      }

      delete read;
      delete compute;
      delete write;

      while(!text_blocks.empty()) {
        delete text_blocks.back();
        text_blocks.pop_back();
      }
    }

    void merge(wavelet_structure_external& result) {
      auto& bvs = wavelet_structure_external_writer::bvs(result);





      DDE2_VERBOSE << "STARTING MERGER. TEXT_SIZE " << size_
                   << " LEVELS " << levels_ << " BV_SIZE " << bvs.size() << "\n";


      const uint64_t number_of_intervals = (1ULL << levels_) - 1;
      std::vector<uint64_t> global_flat_hist(number_of_intervals, 0);
      std::vector<uint64_t> global_flat_borders(number_of_intervals, 0);

      global_flat_hist[0] = size_;
      #pragma omp parallel for
      for (uint64_t i = 1; i < number_of_intervals; ++i) {
        for (uint64_t b = 0; b < number_of_blocks_; ++b) {
          global_flat_hist[i] += block_hists_[b][i];
        }
      }

      global_flat_borders[0] = 0;
//      DDE2_VERBOSE << "0   0   " << size_ << "\n";
      for (uint64_t l = 1; l < levels_; ++l) {
        const uint64_t level_hist_size = 1ULL << l;
        const uint64_t start_idx = level_hist_size - 1;
        const uint64_t end_idx = start_idx + level_hist_size;
        global_flat_borders[start_idx] = global_flat_borders[start_idx - 1] +
                                         global_flat_hist[start_idx - 1];
        global_flat_borders[start_idx] = mul64(div64(global_flat_borders[start_idx] + 63));
//        DDE2_VERBOSE << start_idx << "   " << global_flat_borders[start_idx] << "   " << global_flat_hist[start_idx] << "\n";
        for (uint64_t i = start_idx + 1; i < end_idx; ++i) {
          global_flat_borders[i] = global_flat_borders[i - 1] +
                                   global_flat_hist[i - 1];
//          DDE2_VERBOSE << i << "   " << global_flat_borders[i] << "   " << global_flat_hist[i] << "\n";
        }
      }

      result_reader_type reader(unmerged_result);
      DDE2_VERBOSE << "Initializing buffers. ";
      em_merger merger(levels_, bytes_memory_);
//      em_merger merger(levels_, 1);
      DDE2_VERBOSE << "Total words memory: " << merger.remaining_capacity_in_words() << ".\n";

      for (uint64_t b = 0; b < number_of_blocks_; ++b) {
        const auto& cur_hist = block_hists_[b];
        for (uint64_t i = 0; i < number_of_intervals; ++i) {
          const auto bits = cur_hist[i];
          if (bits > 0) {
            uint64_t words = div64(bits + 63);
            while (PWM_UNLIKELY(words > merger.remaining_capacity_in_words())) {
              const auto largest_interval_idx = merger.largest_interval();
              auto &largest_interval = merger.interval(largest_interval_idx);
              auto border = global_flat_borders[largest_interval_idx];
              global_flat_borders[largest_interval_idx] += largest_interval.size_in_bits();
              result_writer_type writer(bvs.begin() + div64(border));
              // write largest interval to em
              DDE2_VERBOSE << "Write interval " << largest_interval_idx
                           << " of bit length " << largest_interval.size_in_bits()
                           << " to em bit idx " << border << ". Shift = " << mod64(border) << ".\n";
              largest_interval.to_em(writer, mod64(border));
            }

            DDE2_VERBOSE << "Write " << bits << " bits to buffer interval " << i << ".\n";
            // write to buffer
            merger.interval(i).to_memory(reader, bits);
          }
        }
      }

      DDE2_VERBOSE << "Finalize merger. Longest interval: " << merger.largest_interval() << " ("
                   << merger.interval(merger.largest_interval()).size_in_bits() << ").\n";
      for (uint64_t i = 0; i < number_of_intervals; ++i) {
        if (merger.interval(i).size_in_bits() > 0) {
          auto border = global_flat_borders[i];
          result_writer_type writer(bvs.begin() + div64(border));
          DDE2_VERBOSE << "Write interval " << i
                       << " of bit length " << merger.interval(i).size_in_bits()
                       << " to em bit idx " << border << ". Shift = " << mod64(border) << ".\n";
          merger.interval(i).to_em(writer, mod64(border));
        }
      }
    }


    dd_ctx(const InputType &text, const uint64_t size, const uint64_t levels)
        : omp_size_(get_omp_size()),
          text_(text), size_(size), levels_(levels),
          input_bytes_per_block_((bytes_memory_ / omp_size_) >> 1),
          input_chars_per_block_(input_bytes_per_block_ / sizeof(char_type)),
          number_of_blocks_((size + input_chars_per_block_ - 1) / input_chars_per_block_),
          input_chars_last_block_(size_ - input_chars_per_block_ * (number_of_blocks_ - 1)) {

      block_hists_.reserve(number_of_blocks_);
      unmerged_result.reserve((size_ * sizeof(char_type)) / 7);

      DDE2_VERBOSE
        << "Created context for wx_dd_fe "
        << "[blocks: " << number_of_blocks_ << ", "
        << "block_chars: " << input_chars_per_block_ << ", "
        << "last_block_chars: " << input_chars_last_block_  << ", "
        << "omp_size: " << omp_size_
        << "].\n";
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

//    wavelet_structure_external_writer::bvs(result).clear();

    using ctx_t = dd_ctx<InputType>;
    ctx_t ctx(text, size, levels);

    ctx.run();
    stats.phase("merge");
    ctx.merge(result);

//    print_structure(std::cout, result.getInternalStructure(), true);
//
//    auto result2 = wx_dd_pc_fe<typename InputType::value_type, true, bytes_memory_ / 8>::compute(text, size, levels, stats);
//    print_structure(std::cout, result2.getInternalStructure(), true);


    return result;
  }
}; // class wx_ps

/******************************************************************************/
