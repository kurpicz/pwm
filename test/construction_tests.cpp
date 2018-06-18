/*******************************************************************************
 * test/construction_test.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <gtest/gtest.h>

#include "test/util.hpp"

#include "benchmark/algorithm.hpp"
#include "huffman/huff_decode.hpp"
#include "util/alphabet_util.hpp"
#include "util/common.hpp"
#include "util/debug.hpp"
#include "util/file_stream.hpp"
#include "util/file_util.hpp"

template <typename list_type>
void construction_smoketest(list_type& algo_list) {
  for (const auto& a : algo_list) {
    if (a->word_width() == 1 && !a->is_huffman_shaped()) {
      a->print_info();
      test::roundtrip_batch([&](const std::string& s){
        auto vec = std::vector<uint8_t>(s.begin(), s.end());
        uint64_t levels = levels_for_max_char(no_reduction_alphabet(vec));
        auto bvz = a->compute_bitvector(&vec, vec.size() , levels);
        if (a->is_tree()) {
          auto decoded_s = decode_wt(bvz.bvs(), vec.size());
          ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
        } else {
          auto decoded_s = decode_wm(bvz.bvs(), bvz.zeros(), vec.size());
          ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
        }
      });
    }
  }
}

TEST(wavelet_construction, smoketest) {
  auto& algo_list1 = algorithm_list<true>::get_algorithm_list();
  auto& algo_list2 = algorithm_list<false>::get_algorithm_list();
  
  construction_smoketest(algo_list1);
  construction_smoketest(algo_list2);
}

//~ TEST(semi_external_wavelet_construction, smoketest) {
  //~ auto& algo_list = algorithm_list::get_algorithm_list();
  //~ for (const auto& a : algo_list) {
    //~ if (a->word_width() == 1 && !a->is_huffman_shaped()) {
      //~ a->print_info();
      //~ test::roundtrip_batch([&](const std::string& s){
        //~ auto vec = std::vector<uint8_t>(s.begin(), s.end());
        //~ const std::string file_name = "test_string";
        //~ uint64_t levels = levels_for_max_char(no_reduction_alphabet(vec));
        //~ vector_to_file(vec, file_name);
        //~ ifile_stream<uint8_t> ifs(file_name);
        //~ auto bvz = a->compute_bitvector_semi_external(
          //~ &ifs, vec.size() , levels);
        //~ ASSERT_TRUE(remove(file_name.c_str()) == 0) <<
          //~ "Could not remove file, something went wrong";
        //~ if (a->is_tree()) {
          //~ auto decoded_s = decode_wt(bvz.bvs(), vec.size());
          //~ ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
        //~ } else {
          //~ auto decoded_s = decode_wm(bvz.bvs(), bvz.zeros(), vec.size());
          //~ ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
        //~ }
      //~ });
    //~ }
  //~ }
//~ }

template <typename list_type>
void huffman_shaped_wavelet_construction_smoketest(list_type& algo_list) {
  for (const auto& a : algo_list) {
    if (a->word_width() == 1 && a->is_huffman_shaped()) {
      a->print_info();
      test::roundtrip_batch([&](const std::string& s){
        auto vec = std::vector<uint8_t>(s.begin(), s.end());
        uint64_t levels = levels_for_max_char(no_reduction_alphabet(vec));
        auto bvz = a->compute_bitvector(&vec, vec.size() , levels);
        if (a->is_tree()) {
          auto decoded_s = decode_wt_huff(bvz.bvs(), vec.size());
          ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
        } else {
          auto decoded_s = decode_wm_huff(bvz.bvs(), bvz.zeros(), vec.size());
          ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
        }
      });
    }
  }
}

TEST(huffman_shaped_wavelet_construction, smoketest) {
  auto& algo_list1 = algorithm_list<true>::get_algorithm_list();
  auto& algo_list2 = algorithm_list<false>::get_algorithm_list();
  
  huffman_shaped_wavelet_construction_smoketest(algo_list1);
  huffman_shaped_wavelet_construction_smoketest(algo_list2);
}

/******************************************************************************/
