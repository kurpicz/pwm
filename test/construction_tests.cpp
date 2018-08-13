/*******************************************************************************
 * test/construction_test.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 * 
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <bitset>

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
void internal_construction_smoketest(list_type& algo_list) {
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

TEST(internal_wavelet_construction, smoketest) {
  auto& algo_list = 
    algorithm_list<memory_mode::internal>::get_algorithm_list();
  internal_construction_smoketest(algo_list);
}


template <typename list_type>
void external_input_construction_smoketest(list_type& algo_list) {
  for (const auto& a : algo_list) {
    if (a->word_width() == 1 && !a->is_huffman_shaped()) {
      a->print_info();
      test::roundtrip_batch([&](const std::string& s){
        auto vec = std::vector<uint8_t>(s.begin(), s.end());
        uint64_t levels = levels_for_max_char(no_reduction_alphabet(vec));
        
        stxxlvector<type_for_bytes<1>::type> * vec_external = 
          new stxxlvector<type_for_bytes<1>::type>();
        for(const auto symbol : vec) 
          (*vec_external).push_back(symbol);
        
        auto bvz = a->compute_bitvector(vec_external, vec.size() , levels);
        delete vec_external;        
        
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

TEST(external_input_wavelet_construction, smoketest) {
  auto& algo_list = 
    algorithm_list<memory_mode::external_input>::get_algorithm_list();
  external_input_construction_smoketest(algo_list);
}

template <typename list_type>
void external_output_construction_smoketest(list_type& algo_list) {
  for (const auto& a : algo_list) {
    if (a->word_width() == 1 && !a->is_huffman_shaped()) {
      a->print_info();
      test::roundtrip_batch([&](const std::string& s){
        auto vec = std::vector<uint8_t>(s.begin(), s.end());
        uint64_t levels = levels_for_max_char(no_reduction_alphabet(vec));
        auto bvz = a->compute_bitvector(&vec, vec.size(), levels);
        
        if(vec.size() == 0) {
          ASSERT_EQ(bvz.levels(), uint64_t(0)) << "Failure (Algorithm): " << a->name();
        } else {
          internal_bit_vectors bvz_internal(levels, vec.size());
          for(uint64_t level = 0; level < levels; ++level) {
            for(uint64_t entry = 0; entry < (vec.size() + 63) / 64; ++entry) {
              bvz_internal[level][entry] = bvz[level][entry];
            }
          }
          if (a->is_tree()) {
            auto decoded_s = decode_wt(bvz_internal, vec.size());
            ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
          } else {
            if(s.size() > 1 || true) {
              auto decoded_s = decode_wm(bvz_internal, bvz.zeros(), vec.size());
              ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
            }
          }
        }
      });
    }
  }
}

TEST(external_output_wavelet_construction, smoketest) {
  auto& algo_list = 
    algorithm_list<memory_mode::external_output>::get_algorithm_list();
  external_output_construction_smoketest(algo_list);
}

template <typename list_type>
void external_construction_smoketest(list_type& algo_list) {
  for (const auto& a : algo_list) {
    if (a->word_width() == 1 && !a->is_huffman_shaped()) {
      a->print_info();
      test::roundtrip_batch([&](const std::string& s){
        auto vec = std::vector<uint8_t>(s.begin(), s.end());
        uint64_t levels = levels_for_max_char(no_reduction_alphabet(vec));
        
        stxxlvector<type_for_bytes<1>::type> * vec_external = 
          new stxxlvector<type_for_bytes<1>::type>();
        for(const auto symbol : vec) 
          (*vec_external).push_back(symbol);
        
        auto bvz = a->compute_bitvector(vec_external, vec.size() , levels);
        
        if(vec.size() == 0) {
          ASSERT_EQ(bvz.levels(), uint64_t(0)) << "Failure (Algorithm): " << a->name();
        } else {
          internal_bit_vectors bvz_internal(levels, vec.size());
          for(uint64_t level = 0; level < levels; ++level) {
            for(uint64_t entry = 0; entry < (vec.size() + 63) / 64; ++entry) {
              bvz_internal[level][entry] = bvz[level][entry];
            }
          }
          if (a->is_tree()) {
            auto decoded_s = decode_wt(bvz_internal, vec.size());
            ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
          } else {
            if(s.size() > 1 || true) {
              auto decoded_s = decode_wm(bvz_internal, bvz.zeros(), vec.size());
              ASSERT_EQ(s, decoded_s) << "Failure (Algorithm): " << a->name();
            }
          }
        }
      });
    }
  }
}

TEST(external_wavelet_construction, smoketest) {
  auto& algo_list = 
    algorithm_list<memory_mode::external>::get_algorithm_list();
  external_construction_smoketest(algo_list);
}


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
  auto& algo_list1 = algorithm_list<memory_mode::internal>::get_algorithm_list();
  huffman_shaped_wavelet_construction_smoketest(algo_list1);
}

/******************************************************************************/
