/*******************************************************************************
 * src/benchmark.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <tclap/CmdLine.h>
#include <vector>

#include <stxxl/vector>
#include <stxxl/bits/io/linuxaio_file.h>

#include "benchmark/algorithm.hpp"
#include "util/alphabet_util.hpp"
#include "util/file_util.hpp"
#include "util/stxxl_helper.hpp"
#include "util/memory_types.hpp"

#ifdef MALLOC_COUNT
#include "benchmark/malloc_count.h"
#endif // MALLOC_COUNT

auto filter_parallel(bool only_parallel, bool is_parallel) {
  return (!only_parallel || is_parallel);
}

auto filter_sequential(bool sequential, bool is_parallel) {
  return (!sequential || !is_parallel);
}

auto filter_wavelet_type(bool is_tree, bool no_trees, bool no_matrices) {
  return (is_tree ? !no_trees : !no_matrices);
}

static TCLAP::SwitchArg list_all_algorithms("l", "list",
    "Print the name and description of all registered algorithms", false);
static TCLAP::MultiArg<std::string> file_path_arg("f", "file",
    "Path to the text file.", false, "string");
static TCLAP::ValueArg<std::string> filter_arg("n", "name",
    "Runs all algorithms that contain the <name> in their name", false, "",
    "string");
static TCLAP::ValueArg<int32_t> word_width_arg("b", "byte",
    "Bytes per char in the input text.", false, 1, "uint8_t");
static TCLAP::ValueArg<int32_t> nr_runs_arg("r", "runs",
    "Number of repetitions of the construction algorithm.",
    false, 5, "int32_t");
static TCLAP::SwitchArg run_only_parallel_arg("p", "parallel",
    "Run only parallel construction algorithms.", false);
static TCLAP::SwitchArg run_only_sequential_arg("s", "sequential",
    "Run only sequential construction algorithms.", false);
static TCLAP::SwitchArg no_trees_arg("m", "no_trees",
    "Skip all wavelet trees construction algorithms.", false);
static TCLAP::SwitchArg no_matrices_arg("t", "no_matrices",
    "Skip all wavelet matrices construction algorithms.", false);
static TCLAP::SwitchArg memory_arg("", "memory",
    "Compute peak memory during construction.", false);
static TCLAP::SwitchArg external_input_arg("i", "external_input",
    "Run only algorithms that stream the input instead of keeping it in RAM.", false);
static TCLAP::SwitchArg external_output_arg("o", "external_output",
    "Run only algorithms that stream the output instead of keeping it in RAM.", false);
static TCLAP::SwitchArg external_both_arg("e", "external",
    "Run only algorithms that use external memory", false);

template <memory_mode mem_mode>
int32_t run();

int32_t main(int32_t argc, char const* argv[]) {
  TCLAP::CmdLine cmd("Benchmark for wavelet tree (and matrix) construction",
    ' ', "0.2");
  cmd.add(list_all_algorithms);
  cmd.add(file_path_arg);
  cmd.add(filter_arg);
  cmd.add(word_width_arg);
  cmd.add(nr_runs_arg);
  cmd.add(run_only_parallel_arg);
  cmd.add(run_only_sequential_arg);
  cmd.add(no_trees_arg);
  cmd.add(no_matrices_arg);
  cmd.add(memory_arg);
  cmd.add(external_input_arg);
  cmd.add(external_output_arg);
  cmd.add(external_both_arg);
  cmd.parse( argc, argv );

  if(external_both_arg.getValue()) 
    return run<memory_mode::external>();
  else if(external_input_arg.getValue() && external_output_arg.getValue()) 
    return run<memory_mode::external>();
  else if(external_input_arg.getValue()) 
    return run<memory_mode::external_input>();
  else if(external_output_arg.getValue()) 
    return run<memory_mode::external_output>();
  else return run<memory_mode::internal>();
}

template <memory_mode mem_mode>
int32_t run() {

  if (list_all_algorithms.getValue()) {
    auto& algo_list_int = 
      algorithm_list<memory_mode::internal>::get_algorithm_list();
    auto& algo_list_ext_in = 
      algorithm_list<memory_mode::external_input>::get_algorithm_list();
    auto& algo_list_ext_out = 
      algorithm_list<memory_mode::external_output>::get_algorithm_list();
    auto& algo_list_ext_both = 
      algorithm_list<memory_mode::external>::get_algorithm_list();
    for (const auto& a : algo_list_int) a->print_info();
    for (const auto& a : algo_list_ext_in) a->print_info();
    for (const auto& a : algo_list_ext_out) a->print_info();
    for (const auto& a : algo_list_ext_both) a->print_info();
    return 0;
  }
  
  auto& algo_list = algorithm_list<mem_mode>::get_algorithm_list();

  const std::vector<std::string> file_paths = file_path_arg.getValue();
  std::string filter = filter_arg.getValue();
  const int32_t word_width = word_width_arg.getValue();
  const int32_t nr_runs = nr_runs_arg.getValue();
  const bool run_only_parallel = run_only_parallel_arg.getValue();
  const bool run_only_sequential = run_only_sequential_arg.getValue();
  const bool no_trees = no_trees_arg.getValue();
  const bool no_matrices = no_matrices_arg.getValue();
  const bool memory = memory_arg.getValue();

  constexpr bool ext_input = 
    mem_mode == memory_mode::external || 
    mem_mode == memory_mode::external_input;
  constexpr bool ext_output = 
    mem_mode == memory_mode::external || 
    mem_mode == memory_mode::external_output;
  
  for (const auto& path : file_paths) {
    std::cout << std::endl << "Text: " << path << std::endl;
    void * txt_prt = nullptr;
    
    uint64_t txt_bytes = 0;
    uint64_t all_bytes = 0;
    uint64_t text_size = 0;
    uint64_t max_char = 0;
    uint64_t levels = 0;
    std::string txt_path;
    
    std::vector<uint8_t> text_uint8;
    std::vector<uint16_t> text_uint16;
    std::vector<uint32_t> text_uint32;
    std::vector<uint64_t> text_uint64;
    uint8_t * text_uint8_data = nullptr;
    uint16_t * text_uint16_data = nullptr;
    uint32_t * text_uint32_data = nullptr;
    uint64_t * text_uint64_data = nullptr;

    stxxlvector<type_for_bytes<1>::type> * text_uint8_ext = nullptr;
    stxxlvector<type_for_bytes<2>::type> * text_uint16_ext = nullptr;
    stxxlvector<type_for_bytes<4>::type> * text_uint32_ext = nullptr;
    stxxlvector<type_for_bytes<8>::type> * text_uint64_ext = nullptr;

    malloc_count_reset_peak();
    uint64_t non_text_bytes = malloc_count_current();
    
    if(!ext_input) {
      // INTERNAL MEMORY INPUT
      if (word_width == 1) {
        text_uint8 = file_to_vector<1>(path);
        text_size = text_uint8.size();
        max_char = reduce_alphabet(text_uint8);
        text_uint8_data = text_uint8.data();
        txt_prt = &text_uint8_data;
      } else if (word_width == 2) {
        text_uint16 = file_to_vector<2>(path);
        text_size = text_uint16.size();
        max_char = reduce_alphabet(text_uint16);
        text_uint16_data = text_uint16.data();
        txt_prt = &text_uint16_data;
      } else if (word_width == 4) {
        text_uint32 = file_to_vector<4>(path);
        text_size = text_uint32.size();
        max_char = reduce_alphabet(text_uint32);
        text_uint32_data = text_uint32.data();
        txt_prt = &text_uint32_data;
      } else if (word_width == 8) {
        text_uint64 = file_to_vector<8>(path);
        text_size = text_uint64.size();
        max_char = reduce_alphabet(text_uint64);
        text_uint64_data = text_uint64.data();
        txt_prt = &text_uint64_data;
      } else {
        std::cerr << "You entered an invalid number of bytes per character "
                     "(parameter 'b')." << std::endl;
        return -1;
      }
    } else {
      // EXTERNAL MEMORY INPUT
      stxxl::linuxaio_file stxxl_file(path, stxxl::file::open_mode::RDONLY);
      if (word_width == 1) {
        const stxxlvector<type_for_bytes<1>::type> unreduced_vector(&stxxl_file);
        text_size = unreduced_vector.size();
        text_uint8_ext = new stxxlvector<type_for_bytes<1>::type>();
        max_char = reduce_alphabet<type_for_bytes<1>::type>(unreduced_vector, *text_uint8_ext);
        txt_prt = text_uint8_ext;        
      } else if (word_width == 2) {
        const stxxlvector<type_for_bytes<2>::type> unreduced_vector(&stxxl_file);
        text_size = unreduced_vector.size();
        text_uint16_ext = new stxxlvector<type_for_bytes<2>::type>();
        max_char = reduce_alphabet<type_for_bytes<2>::type>(unreduced_vector, *text_uint16_ext);
        txt_prt = text_uint16_ext;        
      } else if (word_width == 4) {
        const stxxlvector<type_for_bytes<4>::type> unreduced_vector(&stxxl_file);
        text_size = unreduced_vector.size();
        text_uint32_ext = new stxxlvector<type_for_bytes<4>::type>();
        max_char = reduce_alphabet<type_for_bytes<4>::type>(unreduced_vector, *text_uint32_ext);
        txt_prt = text_uint32_ext;        
      } else if (word_width == 8) {
        const stxxlvector<type_for_bytes<8>::type> unreduced_vector(&stxxl_file);
        text_size = unreduced_vector.size();
        text_uint64_ext = new stxxlvector<type_for_bytes<8>::type>();
        max_char = reduce_alphabet<type_for_bytes<8>::type>(unreduced_vector, *text_uint64_ext);
        txt_prt = text_uint64_ext;        
      } else {
        std::cerr << "You entered an invalid number of bytes per character "
                     "(parameter 'b')." << std::endl;
        return -1;
      }
    }
    levels = levels_for_max_char(max_char);
    
    #ifdef MALLOC_COUNT
      all_bytes = malloc_count_current();
      txt_bytes = all_bytes - non_text_bytes;
    #endif // MALLOC_COUNT

    std::cout << "Characters: " << text_size << std::endl;
    std::cout << "Levels:     " << levels << std::endl;
    
    #ifdef MALLOC_COUNT
      std::cout << "Memory peak text:  " << txt_bytes << ", MB: "
                << txt_bytes / (1024 * 1024) << std::endl;
      std::cout << "Memory peak total: " << all_bytes << ", MB: "
                << all_bytes / (1024 * 1024) << std::endl;
    #endif // MALLOC_COUNT

    for (const auto& a : algo_list) {
      if (filter == "" || (a->name().find(filter) != std::string::npos)) {
        if (a->word_width() == word_width) {
          if (filter_parallel(run_only_parallel, a->is_parallel())) {
            if (filter_sequential(run_only_sequential, a->is_parallel())) {
              if (filter_wavelet_type(a->is_tree(), no_trees, no_matrices)) {
                a->print_info();
                if (memory) {
                  #ifdef MALLOC_COUNT
                    malloc_count_reset_peak();
                    a->memory_peak(txt_prt, text_size, levels);
                    std::cout << malloc_count_peak() << ", MB: "
                              << malloc_count_peak() / (1024 * 1024) << std::endl;
                  #else
                    std::cout << "Memory measurement is NOT enabled."
                              << std::endl;
                  #endif // MALLOC_COUNT
                } else {
                  std::cout << a->median_time(
                    txt_prt, text_size, levels, nr_runs) << std::endl;
                }
              }
            }
          }
        }
      }
    }
    
    delete text_uint8_ext; 
    delete text_uint16_ext;
    delete text_uint32_ext;
    delete text_uint64_ext;
  }
  return 0;
}

/******************************************************************************/
