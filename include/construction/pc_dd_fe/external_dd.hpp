/*******************************************************************************
 * include/construction/pc_dd_fe/external_merge.hpp
 *
 * Copyright (C) 2019 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/


#pragma once

#include "util/print.hpp"
#include "util/permutation.hpp"
#include "construction/pc_dd_fe/ctx_partial.hpp"
#include "construction/pc_dd_fe/pc_partial.hpp"
#include "construction/pc_dd_fe/external_merge.hpp"

#ifdef WX_DD_PC_FE_VERBOSE
#define EXT_DD_CTX_VERBOSE WX_DD_PC_FE_VERBOSE
#else
#define EXT_DD_CTX_VERBOSE if constexpr (false) atomic_out()
#endif

template <typename InputType, uint64_t bytes_memory, bool rw_simultaneously>
struct external_dd_ctx {
private:

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

  static uint64_t get_omp_size() {
    int result;
    #pragma omp parallel
    result = omp_get_num_threads();
    return (uint64_t)result;
  }

public:
  external_dd_ctx(const InputType &text, const uint64_t size, const uint64_t levels)
      : omp_size_(get_omp_size()),
        text_(text), size_(size), levels_(levels),
        input_bytes_per_block_(std::max((bytes_memory / omp_size_) >> 1, (uint64_t)8)),
        input_chars_per_block_(input_bytes_per_block_ / sizeof(char_type)),
        number_of_blocks_((size + input_chars_per_block_ - 1) / input_chars_per_block_),
        input_chars_last_block_(size_ - input_chars_per_block_ * (number_of_blocks_ - 1)) {

    block_hists_.reserve(number_of_blocks_);
    unmerged_result.reserve((size_ * sizeof(char_type)) / 7);

    EXT_DD_CTX_VERBOSE
        << "Created context for wx_dd_fe "
        << "[blocks: " << number_of_blocks_ << ", "
        << "block_chars: " << input_chars_per_block_ << ", "
        << "last_block_chars: " << input_chars_last_block_  << ", "
        << "omp_size: " << omp_size_
        << "].\n";
  }

  void dd() {
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

          #pragma omp task \
          priority(number_of_blocks_ - b) \
          depend(in: read[b + omp_size_ - 1], compute[b]) \
          depend(out: read[b + omp_size_])
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

            EXT_DD_CTX_VERBOSE << "Read " << b << " (bs " << block_size << ")\n";
            if constexpr (!rw_simultaneously) omp_unset_lock(&rw_lock);
          }

          #pragma omp task \
          priority(0) \
          depend(in: read[b + omp_size_], write[b]), \
          depend(out: compute[b + omp_size_])
          {
            const uint64_t block_size = get_chars_per_block(b);
            char_type * text_block = read[b + omp_size_];
            EXT_DD_CTX_VERBOSE << "START Compute " << b << " " << text_block[1] << " (bs " << block_size << ")\n";
            ctx_partial * ctx_block = new ctx_partial(block_size, levels_);
            pc_partial(text_block, ctx_block);
            compute[b + omp_size_] = ctx_block;
            EXT_DD_CTX_VERBOSE << "DONE Compute " << b << " " << text_block[1] << " (bs " << block_size << ")\n";
          }

          #pragma omp task \
          priority(number_of_blocks_ - b) \
          depend(in: compute[b + omp_size_], write[b + omp_size_ - 1]) \
          depend(out: write[b + omp_size_])
          {
            if constexpr (!rw_simultaneously) omp_set_lock(&rw_lock);
            EXT_DD_CTX_VERBOSE << "START Write " << b << " (bs " << get_chars_per_block(b) << ")\n";
            ctx_partial * ctx_block = compute[b + omp_size_];
            block_hists_.emplace_back(ctx_block->flat_hist());
            const auto &block_bvs = ctx_block->bv();
            const auto &level_data_sizes = ctx_block->level_data_sizes();
            EXT_DD_CTX_VERBOSE << "Levels " << levels_ << ", ds " << ctx_block->data_size() << "\n";
            for (uint64_t l = 0; l < levels_; ++l) {
              const auto level_data_size = level_data_sizes[l];
              EXT_DD_CTX_VERBOSE << "    Level " << l << ", ds " << level_data_size << "\n";
              for (uint64_t w = 0; w < level_data_size; ++w) {
                writer << block_bvs[l][w];
                EXT_DD_CTX_VERBOSE << std::bitset<64>(block_bvs[l][w]).to_string() << "\n";
              }
            }
            delete ctx_block;
            EXT_DD_CTX_VERBOSE << "DONE Write " << b << " (bs " << get_chars_per_block(b) << ")\n";
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

  template <bool is_tree>
  void merge(wavelet_structure_external& result) {
    auto& bvs = wavelet_structure_external_writer::bvs(result);


    EXT_DD_CTX_VERBOSE << "STARTING MERGER. TEXT_SIZE " << size_
                       << " LEVELS " << levels_ << " BV_SIZE "
                       << bvs.size() << " IS TREE " << is_tree << "\n";

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

    auto rho = rho_dispatch<is_tree>::create(levels_);

    global_flat_borders[0] = 0;

    uint64_t s = 1;
    for (uint64_t l = 1; l < levels_; ++l) {
      const uint64_t level_hist_size = 1ULL << l;
      global_flat_borders[s] = global_flat_borders[s - 1] +
                                   global_flat_hist[s - 1];
      global_flat_borders[s] = mul64(div64(global_flat_borders[s] + 63));

      for (uint64_t i = 1; i < level_hist_size; ++i) {
        global_flat_borders[s + rho(l, i)] =
            global_flat_borders[s + rho(l, i - 1)] +
            global_flat_hist[s + rho(l, i - 1)];
      }
      s += level_hist_size;
    }

    if constexpr (!is_tree) {
      auto& zeros = wavelet_structure_external_writer::zeros(result);
      uint64_t s = 1;
      for (uint64_t l = 1; l <= levels_; ++l) {
        const uint64_t level_hist_size = 1ULL << l;
        for (uint64_t i = 0; i < level_hist_size; i += 2) {
          zeros[l - 1] += global_flat_hist[s + i];
        }
        s += level_hist_size;
      }
    }


    result_reader_type reader(unmerged_result);
    EXT_DD_CTX_VERBOSE << "Initializing buffers. ";
    em_merger merger(levels_, bytes_memory);
//      em_merger merger(levels_, 1);
    EXT_DD_CTX_VERBOSE << "Total words memory: " << merger.remaining_capacity_in_words() << ".\n";

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
            EXT_DD_CTX_VERBOSE << "Write interval " << largest_interval_idx
                               << " of bit length " << largest_interval.size_in_bits()
                               << " to em bit idx " << border << ". Shift = " << mod64(border) << ".\n";
            largest_interval.to_em(writer, mod64(border));
          }

          EXT_DD_CTX_VERBOSE << "Write " << bits << " bits to buffer interval " << i << ".\n";
          // write to buffer
          merger.interval(i).to_memory(reader, bits);
        }
      }
    }

    EXT_DD_CTX_VERBOSE << "Finalize merger. Longest interval: " << merger.largest_interval() << " ("
                       << merger.interval(merger.largest_interval()).size_in_bits() << ").\n";
    for (uint64_t i = 0; i < number_of_intervals; ++i) {
      if (merger.interval(i).size_in_bits() > 0) {
        auto border = global_flat_borders[i];
        result_writer_type writer(bvs.begin() + div64(border));
        EXT_DD_CTX_VERBOSE << "Write interval " << i
                           << " of bit length " << merger.interval(i).size_in_bits()
                           << " to em bit idx " << border << ". Shift = " << mod64(border) << ".\n";
        merger.interval(i).to_em(writer, mod64(border));
      }
    }
  }
};















