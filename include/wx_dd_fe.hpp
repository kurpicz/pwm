/*******************************************************************************
 * include/wx_ps_fe.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>

#include "arrays/memory_types.hpp"
#include "construction/merge_external.hpp"
#include "construction/wavelet_structure.hpp"
#include "construction/wavelet_structure_external.hpp"

#include "wx_base.hpp"
#include "wx_dd_pc.hpp"

template <typename AlphabetType, bool is_tree_>
class wx_dd_fe : public wx_in_out_external<true, true> {
  static constexpr uint64_t bytesPerBlock = 537;
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
    auto& bvs = wavelet_structure_external_writer::bvs(result);
    bvs.clear();

    // typedefs
    using parallel_algo = wx_dd_pc<AlphabetType, is_tree_>;

    using bv_type = typename std::remove_reference<decltype(bvs)>::type;
    using bv_writer_type = bv_type::bufwriter_type;
    using text_reader_type = typename InputType::bufreader_type;

    // calculate block sizes etc
    const uint64_t block_size = bytesPerBlock / sizeof(AlphabetType) / 64 * 64;
    const uint64_t block_count = (size + block_size - 1) / block_size;
    const uint64_t last_block_size = size - ((block_count - 1) * block_size);

    const uint64_t block_level_words = (block_size + 63) / 64;
    const uint64_t block_words = block_level_words * levels;
    const uint64_t last_block_level_words = (last_block_size + 63) / 64;
    const uint64_t last_block_words = last_block_level_words * levels;
    const uint64_t temp_data_size =
        block_words * (block_count - 1) + last_block_words;


    std::cout << std::endl << name.str()
              << ": block_size=" << block_size
              << ", block_count=" << block_count
              << ", last_block_size=" << last_block_size
              << std::endl;

    // readers, writers, buffers
    bv_type temp_vector;
    temp_vector.reserve(temp_data_size);
    bv_writer_type temp_writer(temp_vector);
    bv_writer_type result_writer(bvs);
    text_reader_type text_reader(text);
    AlphabetType block [block_size];

    // initialize block histograms
    std::vector<std::vector<std::vector<uint64_t>>> block_hists;
    block_hists.resize(block_count);
    for(uint64_t b = 0; b < block_count; b++) {
      block_hists[b].resize(levels + 1);
      block_hists[b][0].resize(1, 0);
      for(uint64_t l = 1; l < levels + 1; l++) {
        block_hists[b][l].resize(block_hists[b][l - 1].size() * 2, 0);
      }
    }

    for(uint64_t b = 0; b < block_count; b++) {
      const uint64_t current_block_size =
          (b + 1 < block_count) ? block_size : last_block_size;
      const uint64_t current_block_level_words =
          (b + 1 < block_count) ? block_level_words : last_block_level_words;
      //read block to internal memory and calculate last level histograms
      for(uint64_t i = 0; i < current_block_size; i++) {
        auto symbol = *text_reader;
        block[i] = symbol;
        ++block_hists[b][levels][symbol];
        ++text_reader;
      }
      // calculate remaining histograms
      for(uint64_t l = levels - 1; l > 0; l--) {
        for(uint64_t s = 0; s < block_hists[b][l].size(); s++) {
          block_hists[b][l][s] =
              block_hists[b][l + 1][2 * s] + block_hists[b][l + 1][2 * s + 1];
        }
      }

      auto block_ws = parallel_algo::compute(&(block[0]),
                                             current_block_size,
                                             levels);
      auto& block_bvs = block_ws.bvs();

      for(uint64_t l = 0; l < levels; l++) {
        auto& writer = (l == 0) ? result_writer : temp_writer;
        for(uint64_t w = 0; w < current_block_level_words; w++) {
          writer << block_bvs[l][w];
//          std::cout << block_bvs[l][w] << std::endl;
        }
      }
    }
    temp_writer.finish();

    std::vector<uint64_t> new_block_offsets(block_count, 0);
    for(uint64_t b = 1; b < block_count; b++)
      new_block_offsets[b] = new_block_offsets[b - 1] + block_size;

    external_merger merger(result_writer, temp_vector);

    //MERGE:
    for(uint64_t l = 1; l < levels; l++) {
      auto block_offsets = new_block_offsets;
      for(uint64_t s = 0; s < block_hists[0][l].size(); s++) {
        for (uint64_t b = 0; b < block_count; b++) {
          auto bits_to_copy = block_hists[b][l][s];
          auto block_offset = block_offsets[b];

          // merge

          block_offsets[b] += bits_to_copy;
        }
      }
    }


    result_writer.finish();

    return result;
  }
}; // class wx_ps

/******************************************************************************/
