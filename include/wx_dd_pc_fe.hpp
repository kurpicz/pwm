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
#include "construction/merge_external.hpp"
#include "construction/wavelet_structure.hpp"
#include "construction/wavelet_structure_external.hpp"
#include "util/permutation.hpp"

#include "wx_base.hpp"
#include "wx_dd_pc.hpp"

#define DDE_VERBOSE if constexpr (false) atomic_out()


template <typename AlphabetType, bool is_tree_, uint64_t bytesPerBlock = 1024 * 1024 * 1024>
class wx_dd_pc_fe : public wx_in_out_external<true, true> {
  static constexpr uint64_t given_block_chars = bytesPerBlock / sizeof(AlphabetType);
  static constexpr uint64_t max_block_chars =
      std::max(uint64_t(64), uint64_t(given_block_chars / 64 * 64));

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

  template <typename InputType, bool dd_verbose = false>
  struct dd_ctx {
    using result_type = typename wavelet_structure_external::bv_type;
    using result_writer_type = typename result_type::bufwriter_type;
    using text_reader_type = typename InputType::bufreader_type;
    using parallel_algo = wx_dd_pc<AlphabetType, is_tree_>;
    using sequential_algo = wx_pc<AlphabetType, is_tree_>;

    const uint64_t omp_size;
    const uint64_t levels;

    const uint64_t block_chars;
    const uint64_t block_count;
    const uint64_t last_block_chars;

    const uint64_t block_level_words;
    const uint64_t block_level_bits;
    const uint64_t last_block_level_words;
    const uint64_t last_block_level_bits;
    const uint64_t result_level_words;

    const uint64_t block_words;
    const uint64_t block_bits;
    const uint64_t result_words;

    const std::vector<uint64_t> bit_reverse;

    std::vector<AlphabetType> block_vec1;
    std::vector<AlphabetType> block_vec2;
    AlphabetType * frontIn = block_vec1.data();
    AlphabetType * backIn = block_vec2.data();

    std::vector<std::vector<std::vector<uint64_t>>> block_hists;

    wavelet_structure * frontRes = nullptr;
    wavelet_structure * backRes = nullptr;

    text_reader_type text_reader;

    result_type temp_result;
    result_writer_type temp_result_writer;

    dd_ctx(const InputType &text, const uint64_t size, const uint64_t plevels)
        : omp_size(omp_get_num_threads()),
          levels(plevels),
          block_chars(std::min(max_block_chars, size)),
          block_count((size + block_chars - 1) / block_chars),
          last_block_chars(size - ((block_count - 1) * block_chars)),
          block_level_words((block_chars + 63) / 64),
          block_level_bits(block_level_words * 64),
          last_block_level_words((last_block_chars + 63) / 64),
          last_block_level_bits(last_block_level_words * 64),
          result_level_words(block_level_words * (block_count - 1) + last_block_level_words),
          block_words(block_level_words * levels),
          block_bits(block_words * 64),
          result_words(result_level_words * levels),
          bit_reverse(bit_reverse_permutation(levels - 1)),
          block_vec1(block_chars),
          block_vec2(block_chars),
          frontIn(block_vec1.data()),
          backIn(block_vec2.data()),
          block_hists(block_count),
          text_reader(text),
          temp_result_writer(temp_result) {
      temp_result = stxxl_files::getVectorTemporary<result_type>(1);
      temp_result.reserve(result_words);
      DDE_VERBOSE
          << "Created context for wx_dd_fe "
          << "[blocks: " << block_count << ", "
          << "block_chars: " << block_chars << ", "
          << "last_block_chars: " << last_block_chars
          << "].\n";
    }

    ~dd_ctx() {
      delete frontRes;
      delete backRes;
    }

    inline void swapInputs() {
      std::swap(frontIn, backIn);
    }

    inline void swapResults() {
      std::swap(frontRes, backRes);
    }

    inline void loadBackInput(uint64_t b) {
      DDE_VERBOSE << "(load " << b << " start) ";

      block_hists[b].resize(levels + 1);
      block_hists[b][0].resize(1, 0);
      for(uint64_t l = 1; l < levels + 1; l++) {
        block_hists[b][l].resize(block_hists[b][l - 1].size() * 2, 0);
      }

      const uint64_t current_block_chars =
          (b + 1 < block_count) ?
          block_chars : last_block_chars;

      //read block to internal memory
      for(uint64_t i = 0; i < current_block_chars; i++) {
        text_reader >> backIn[i];
      }
      DDE_VERBOSE << "(load " << b << " done) ";
    }

    inline void processFrontInput(uint64_t b) {
      DDE_VERBOSE << "(process " << b << " start) ";
      const uint64_t current_block_chars =
          (b + 1 < block_count) ?
          block_chars : last_block_chars;

      delete frontRes;
      auto& current_block_hist = block_hists[b];

      frontRes = new wavelet_structure(
          (omp_size == 1) ?
              sequential_algo::compute(
                  frontIn,
                  current_block_chars,
                  levels,
                  &current_block_hist[levels])
              :
              parallel_algo::compute(
                  frontIn,
                  current_block_chars,
                  levels,
                  &current_block_hist[levels]));

      if constexpr (is_tree_) {
        // calculate remaining histograms
        for (int64_t l = levels - 1; l >= 0; l--) {
          auto &block_hist_level = current_block_hist[l];
          auto &block_hist_level_plus = current_block_hist[l + 1];
          for (uint64_t s = 0; s < block_hist_level.size(); s++) {
            block_hist_level[s] =
                block_hist_level_plus[2 * s] +
                block_hist_level_plus[2 * s + 1];
          }
        }
      } else {
        {
          auto &block_hist_level = current_block_hist[levels - 1];
          auto &block_hist_level_plus = current_block_hist[levels];
          uint64_t block_hist_level_size = block_hist_level.size();
          for (uint64_t s = 0; s < block_hist_level_size; s++) {
            auto rho_s = bit_reverse[s];
            block_hist_level[rho_s] =
                block_hist_level_plus[2 * s] +
                block_hist_level_plus[2 * s + 1];
          }
        }
        for (int64_t l = levels - 2; l >= 0; l--) {
          auto &block_hist_level = current_block_hist[l];
          auto &block_hist_level_plus = current_block_hist[l + 1];
          uint64_t block_hist_level_size = block_hist_level.size();
          for (uint64_t s = 0; s < block_hist_level_size; s++) {
            block_hist_level[s] =
                block_hist_level_plus[s] +
                block_hist_level_plus[s + block_hist_level_size];
          }
        }
      }
      DDE_VERBOSE << "(process " << b << " done) ";
    }

    inline void saveBackResult(uint64_t b) {
      DDE_VERBOSE << "(save " << b << " start) ";
      const uint64_t current_block_level_words =
          (b + 1 < block_count) ?
          block_level_words : last_block_level_words;

      auto& block_bvs = (*backRes).bvs();

      for(uint64_t l = 0; l < levels; l++) {
        for(uint64_t w = 0; w < current_block_level_words; w++) {
          temp_result_writer << block_bvs[l][w];
        }
      }
      DDE_VERBOSE << "(save " << b << " done) ";
    }

    inline void finish() {
      temp_result_writer.finish();
      DDE_VERBOSE << "\n";
    }


    inline void merge(wavelet_structure_external& result) {
      DDE_VERBOSE << "Merging... ";

      auto& bvs = wavelet_structure_external_writer::bvs(result);
      auto& zeros = wavelet_structure_external_writer::zeros(result);
      auto& hists = wavelet_structure_external_writer::histograms(result);

      external_merger merger(temp_result, bvs);
      for(uint64_t l = 0; l < levels; l++) {
        std::vector<uint64_t> block_offsets(block_count);
        for(uint64_t b = 0; b < block_count - 1; b++)
          block_offsets[b] = b * block_bits + l * block_level_bits;
        block_offsets[block_count - 1] =
            (block_count - 1) * block_bits +
            l * last_block_level_bits;
        for(uint64_t s = 0; s < block_hists[0][l].size(); s++) {
          for (uint64_t b = 0; b < block_count; b++) {
            auto bits_to_copy = block_hists[b][l][s];
            auto block_offset = block_offsets[b];
            merger.write(bits_to_copy, block_offset);
            block_offsets[b] += bits_to_copy;
          }
        }
        merger.finishLevel();
      }
      merger.finish();

      for (uint64_t b = 0; b < block_count; ++b) {
        if constexpr (is_tree_) hists[0][0] += block_hists[b][0][0];
        for (uint64_t l = 1; l < levels; ++l) {
          auto &block_hist_level = block_hists[b][l];
          const uint64_t hist_level_size = block_hist_level.size();
          if constexpr (is_tree_) {
            for (uint64_t s = 0; s < hist_level_size; s++) {
              hists[l][s] += block_hist_level[s];
//              if(s % 2 == 0) zeros[l - 1] += block_hist_level[s];
            }
          } else {
            for (uint64_t s = 0; s < hist_level_size / 2; s++) {
//              global_hist_level[s] += block_hist_level[s];
              zeros[l - 1] += block_hist_level[s];
            }
//            for (uint64_t s = hist_level_size / 2; s < hist_level_size; s++) {
//              global_hist_level[s] += block_hist_level[s];
//            }
          }
        }
      }
      DDE_VERBOSE << "Done.\n";
    }
  };

public:
  static constexpr bool is_parallel = true;
  static constexpr bool is_tree = is_tree_;
  static constexpr uint8_t word_width = sizeof(AlphabetType);
  static constexpr bool is_huffman_shaped = false;

  template <typename InputType>
  static wavelet_structure_external
  compute(const InputType& text,
          const uint64_t size,
          const uint64_t levels) {
    static_assert(sizeof(typename InputType::value_type) == sizeof(AlphabetType));

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

    if(ctx.block_count < 2) {
      for(uint64_t b = 0; b < ctx.block_count; b++) {
        ctx.loadBackInput(b);
        ctx.swapInputs();
        ctx.processFrontInput(b);
        ctx.swapResults();
        ctx.saveBackResult(b);
      }
    } else {
      // initialize:
      //   load first input block
      //   load second input block (background task)
      ctx.loadBackInput(0);
      ctx.swapInputs();
      std::thread worker_loadBack([&ctx](){ctx.loadBackInput(1);});
      // no partial result to write to disk yet
      std::thread worker_saveBack([](){});

      for(uint64_t b = 0; b < ctx.block_count - 2; b++) {
        // process current block
        ctx.processFrontInput(b);

        // wait until next input block is loaded from disk
        worker_loadBack.join();
        ctx.swapInputs();
        // load second next input block from disk (background task)
        worker_loadBack =
            std::thread([&ctx, b] () {ctx.loadBackInput(b + 2);});

        // wait until previous result is written to disk
        worker_saveBack.join();
        ctx.swapResults();
        // save current result to disk (background task)
        worker_saveBack =
            std::thread([&ctx, b] () {ctx.saveBackResult(b);});
      }

      ctx.processFrontInput(ctx.block_count - 2);
      worker_loadBack.join();
      ctx.swapInputs();
      worker_saveBack.join();
      ctx.swapResults();
      // save second last result to disk (background task)
      worker_saveBack =
          std::thread([&ctx](){ctx.saveBackResult(ctx.block_count - 2);});
      ctx.processFrontInput(ctx.block_count - 1);
      worker_saveBack.join();
      ctx.swapResults();
      ctx.saveBackResult(ctx.block_count - 1);
    }
    ctx.finish();

    //MERGE:
    ctx.merge(result);

    return result;
  }
}; // class wx_ps

/******************************************************************************/
