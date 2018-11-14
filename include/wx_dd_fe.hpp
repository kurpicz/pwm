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

#include "arrays/memory_types.hpp"
#include "construction/merge_external.hpp"
#include "construction/wavelet_structure.hpp"
#include "construction/wavelet_structure_external.hpp"

#include "wx_base.hpp"
#include "wx_dd_pc.hpp"

template <typename AlphabetType, bool is_tree_, uint64_t bytesPerBlock = 100>
class wx_dd_fe : public wx_in_out_external<true, true> {
  static constexpr uint64_t maxBlockSize = bytesPerBlock / sizeof(AlphabetType) / 64 * 64;

  template <typename InputType>
  struct dd_ctx {
    using result_type = typename wavelet_structure_external::bv_type;
    using result_writer_type = typename result_type::bufwriter_type;
    using text_reader_type = typename InputType::bufreader_type;
    using parallel_algo = wx_dd_pc<AlphabetType, is_tree_>;

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

    std::vector<AlphabetType> block_vec1;
    std::vector<AlphabetType> block_vec2;
    AlphabetType * frontIn = block_vec1.data();
    AlphabetType * backIn = block_vec2.data();

    wavelet_structure * frontRes = nullptr;
    wavelet_structure * backRes = nullptr;

    text_reader_type text_reader;

    result_type result;
    result_writer_type result_writer;

    std::vector<std::vector<std::vector<uint64_t>>> block_hists;

    uint64_t next_block = 0;
    uint64_t current_block = 0;

    dd_ctx(const InputType &text, const uint64_t size, const uint64_t plevels)
        : levels(plevels),
          block_chars(std::min(maxBlockSize, size)),
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
          block_vec1(block_chars),
          block_vec2(block_chars),
          frontIn(block_vec1.data()),
          backIn(block_vec2.data()),
          text_reader(text),
          result_writer(result) {

      block_hists.resize(block_count);
      for(uint64_t b = 0; b < block_count; b++) {
        block_hists[b].resize(levels + 1);
        block_hists[b][0].resize(1, 0);
        for(uint64_t l = 1; l < levels + 1; l++) {
          block_hists[b][l].resize(block_hists[b][l - 1].size() * 2, 0);
        }
      }
      result.reserve(result_words);
    }

    inline void swapInputs() {
      std::swap(frontIn, backIn);
    }

    inline void swapResults() {
      std::swap(frontRes, backRes);
    }

    inline void loadBackInput() {
      current_block = next_block++;
      const uint64_t current_block_chars =
          (next_block < block_count) ?
          block_chars : last_block_chars;

      //read block to internal memory and calculate last level histograms
      for(uint64_t i = 0; i < current_block_chars; i++) {
        auto symbol = *text_reader;
        backIn[i] = symbol;
        ++block_hists[current_block][levels][symbol];
        ++text_reader;
      }
      // calculate remaining histograms
      for(int64_t l = levels - 1; l >= 0; l--) {
        for(uint64_t s = 0; s < block_hists[current_block][l].size(); s++) {
          block_hists[current_block][l][s] =
              block_hists[current_block][l + 1][2 * s] +
              block_hists[current_block][l + 1][2 * s + 1];
        }
      }
    }

    inline void processFrontInput() {
      const uint64_t current_block_chars =
          (next_block < block_count) ?
          block_chars : last_block_chars;

      delete frontRes;

      frontRes = new wavelet_structure(
          parallel_algo::compute(
          frontIn,
          current_block_chars,
          levels));
    }

    inline void saveBackResult() {
      const uint64_t current_block_level_words =
          (next_block < block_count) ?
          block_level_words : last_block_level_words;

      auto& block_bvs = (*backRes).bvs();

      for(uint64_t l = 0; l < levels; l++) {
        for(uint64_t w = 0; w < current_block_level_words; w++) {
          result_writer << block_bvs[l][w];
        }
      }
    }

    inline void finish() {
      assert(next_block == block_count);
      result_writer.finish();
    }


    inline void merge(result_type& bvs) {
      external_merger merger(&result, bvs);
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
            histograms().zeros().
            construct(size, levels, name.str(), 0);
    if(size == 0) return result;

    auto& bvs = wavelet_structure_external_writer::bvs(result);
    bvs.clear();

    const dd_ctx ctx(text, size, levels);

    std::cout << " --INITIAL SCAN-- " << std::flush;
    for(uint64_t b = 0; b < ctx.block_count; b++) {
      ctx.loadBackInput();
      ctx.swapInputs();
      ctx.processFrontInput();
      ctx.swapResults();
      ctx.saveBackResult();
    }
    ctx.finish();

    std::cout << " --MERGE-- " << std::flush;


    //MERGE:
    ctx.merge(bvs);

    return result;
  }
}; // class wx_ps

/******************************************************************************/
