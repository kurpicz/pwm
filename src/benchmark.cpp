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

int32_t main(int32_t argc, char const* argv[]) {
  TCLAP::CmdLine cmd("Benchmark for wavelet tree (and matrix) construction",
    ' ', "0.2");

  TCLAP::SwitchArg list_all_algorithms("l", "list",
    "Print the name and description of all registered algorithms", false);
  cmd.add(list_all_algorithms);
  TCLAP::MultiArg<std::string> file_path_arg("f", "file",
    "Path to the text file.", false, "string");
  cmd.add(file_path_arg);
  TCLAP::ValueArg<std::string> filter_arg("n", "name",
    "Runs all algorithms that contain the <name> in their name", false, "",
    "string");
  cmd.add(filter_arg);
  TCLAP::ValueArg<int32_t> word_width_arg("b", "byte",
    "Bytes per char in the input text.", false, 1, "uint8_t");
  cmd.add(word_width_arg);
  TCLAP::ValueArg<int32_t> nr_runs_arg("r", "runs",
    "Number of repetitions of the construction algorithm.",
    false, 5, "int32_t");
  cmd.add(nr_runs_arg);
  TCLAP::SwitchArg run_only_parallel_arg("p", "parallel",
    "Run only parallel construction algorithms.", false);
  cmd.add(run_only_parallel_arg);
  TCLAP::SwitchArg run_only_sequential_arg("s", "sequential",
    "Run only sequential construction algorithms.", false);
  cmd.add(run_only_sequential_arg);
  TCLAP::SwitchArg no_trees_arg("m", "no_trees",
    "Skip all wavelet trees construction algorithms.", false);
  cmd.add(no_trees_arg);
  TCLAP::SwitchArg no_matrices_arg("t", "no_matrices",
    "Skip all wavelet matrices construction algorithms.", false);
  cmd.add(no_matrices_arg);
  TCLAP::SwitchArg memory_arg("", "memory",
    "Compute peak memory during construction.", false);
  cmd.add(memory_arg);
  TCLAP::SwitchArg external_input_arg("ie", "external_input",
    "Run only algorithms that stream the input instead of keeping it in RAM.", false);
  cmd.add(external_input_arg);
  TCLAP::SwitchArg external_output_arg("oe", "external_output",
    "Run only algorithms that stream the output instead of keeping it in RAM.", false);
  cmd.add(external_output_arg);
  TCLAP::SwitchArg external_both_arg("e", "external",
    "Run only algorithms that use external memory", false);
  cmd.add(external_both_arg);
  
  cmd.parse( argc, argv );

  auto& algo_list_int = 
    algorithm_list<memory_mode::internal>::get_algorithm_list();
  auto& algo_list_ext_in = 
    algorithm_list<memory_mode::external_input>::get_algorithm_list();
  auto& algo_list_ext_out = 
    algorithm_list<memory_mode::external_output>::get_algorithm_list();
  auto& algo_list_ext_both = 
    algorithm_list<memory_mode::external>::get_algorithm_list();
    
  if (list_all_algorithms.getValue()) {
    for (const auto& a : algo_list_int) a->print_info();
    for (const auto& a : algo_list_ext_in) a->print_info();
    for (const auto& a : algo_list_ext_out) a->print_info();
    for (const auto& a : algo_list_ext_both) a->print_info();
    return 0;
  }

  const std::vector<std::string> file_paths = file_path_arg.getValue();
  std::string filter = filter_arg.getValue();
  const int32_t word_width = word_width_arg.getValue();
  const int32_t nr_runs = nr_runs_arg.getValue();
  const bool run_only_parallel = run_only_parallel_arg.getValue();
  const bool run_only_sequential = run_only_sequential_arg.getValue();
  const bool no_trees = no_trees_arg.getValue();
  const bool no_matrices = no_matrices_arg.getValue();
  const bool memory = memory_arg.getValue();
  
  memory_mode mem_mode = memory_mode::internal;  
  if(external_input_arg.getValue()) 
    mem_mode = memory_mode::external_input;
  if(external_output_arg.getValue()) 
    mem_mode = memory_mode::external_output;
  if(external_input_arg.getValue() && external_output_arg.getValue()) 
    mem_mode = memory_mode::external;
  if(external_both_arg.getValue()) 
    mem_mode = memory_mode::external;

  for (const auto& path : file_paths) {
    std::cout << std::endl << "Text: " << path << std::endl;
    void* txt_prt = nullptr;
    void* txt_prt_ext = nullptr;
    uint64_t txt_bytes = 0;
    uint64_t txt_bytes_ext = 0;
    uint64_t text_size = 0;
    uint64_t max_char = 0;
    uint64_t levels = 0;
    std::string txt_path;
    std::vector<uint8_t> text_uint8;
    std::vector<uint16_t> text_uint16;
    std::vector<uint32_t> text_uint32;
    std::vector<uint64_t> text_uint64;
    uint8_t * text_uint8_data;
    uint16_t * text_uint16_data;
    uint32_t * text_uint32_data;
    uint64_t * text_uint64_data;
#ifdef MALLOC_COUNT
    malloc_count_reset_peak();
#endif
    if (false) {
      if (word_width == 1) {
        auto reduced_result = reduce_alphabet_stream<uint8_t>(path);
        txt_path = reduced_result.file_name;
        text_size = reduced_result.file_size;
        max_char = reduced_result.max_char;
      } else if (word_width == 2) {
        auto reduced_result = reduce_alphabet_stream<uint16_t>(path);
        txt_path = reduced_result.file_name;
        text_size = reduced_result.file_size;
        max_char = reduced_result.max_char;
      } else if (word_width == 4) {
        auto reduced_result = reduce_alphabet_stream<uint32_t>(path);
        txt_path = reduced_result.file_name;
        text_size = reduced_result.file_size;
        max_char = reduced_result.max_char;
      } else if (word_width == 8) {
        auto reduced_result = reduce_alphabet_stream<uint64_t>(path);
        txt_path = reduced_result.file_name;
        text_size = reduced_result.file_size;
        max_char = reduced_result.max_char;
      } else {
        std::cerr << "You entered an invalid number of bytes per character "
                     "(parameter 'b')." << std::endl;
        return -1;
      }
    } else {
      
      uint32_t bytes = malloc_count_current();
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
        txt_prt = &text_uint8_data;
      } else if (word_width == 4) {
        text_uint32 = file_to_vector<4>(path);
        text_size = text_uint32.size();
        max_char = reduce_alphabet(text_uint32);
        text_uint32_data = text_uint32.data();
        txt_prt = &text_uint8_data;
      } else if (word_width == 8) {
        text_uint64 = file_to_vector<8>(path);
        text_size = text_uint64.size();
        max_char = reduce_alphabet(text_uint64);
        text_uint64_data = text_uint64.data();
        txt_prt = &text_uint8_data;
      } else {
        std::cerr << "You entered an invalid number of bytes per character "
                     "(parameter 'b')." << std::endl;
        return -1;
      }
      txt_bytes = malloc_count_current() - bytes;
     
      stxxl::linuxaio_file * stxxl_file = 
        new stxxl::linuxaio_file(path, stxxl::file::open_mode::RDONLY);
      
      bytes = malloc_count_current();
      if (word_width == 1) {
        stxxlvector<type_for_bytes<1>::type> * vec = new stxxlvector<type_for_bytes<1>::type>(); 
        for(const auto symbol : text_uint8) { (*vec).push_back(symbol); }
        txt_prt_ext = vec;        
      } else if (word_width == 2) {
        txt_prt_ext = new stxxlvector<type_for_bytes<2>::type>(stxxl_file);
      } else if (word_width == 4) {
        txt_prt_ext = new stxxlvector<type_for_bytes<4>::type>(stxxl_file);
      } else if (word_width == 8) {
        txt_prt_ext = new stxxlvector<type_for_bytes<8>::type>(stxxl_file);
      } else {
        std::cerr << "You entered an invalid number of bytes per character "
                     "(parameter 'b')." << std::endl;
        return -1;
      }
      txt_bytes_ext = malloc_count_current() - bytes;
    }
    levels = levels_for_max_char(max_char);
    std::cout << "Characters: " << text_size << std::endl;
#ifdef MALLOC_COUNT
    std::cout << "Memory peak text: " << malloc_count_peak() - txt_bytes_ext << ", MB: "
              << (malloc_count_peak() - txt_bytes_ext) / (1024 * 1024) << std::endl;
#endif // MALLOC_COUNT
    for (const auto& a : algo_list_int) {
      if (filter == "" || (a->name().find(filter) != std::string::npos)) {
        if (a->word_width() == word_width) {
          if (filter_parallel(run_only_parallel, a->is_parallel())) {
            if (filter_sequential(run_only_sequential, a->is_parallel())) {
              if (filter_wavelet_type(a->is_tree(), no_trees, no_matrices)) {
                a->print_info();
                if (memory) {
#ifdef MALLOC_COUNT
                  malloc_count_reset_peak();
                  if(a->is_input_external()) {
                    a->memory_peak(txt_prt_ext, text_size, levels);
                    std::cout << malloc_count_peak() - txt_bytes << ", MB: "
                            << (malloc_count_peak() - txt_bytes) / (1024 * 1024) << std::endl;
                    std::cout << malloc_count_peak() << ", MB: "
                            << (malloc_count_peak()) / (1024 * 1024) << std::endl;
                  } else {
                    a->memory_peak(txt_prt, text_size, levels);
                    std::cout << malloc_count_peak() - txt_bytes_ext << ", MB: "
                            << (malloc_count_peak() - txt_bytes_ext) / (1024 * 1024) << std::endl;
                  }
#else
                  std::cout << "Memory measurement is NOT enabled."
                            << std::endl;
#endif // MALLOC_COUNT
                } else {
                  if(a->is_input_external())
                    std::cout << a->median_time(
                      txt_prt_ext, text_size, levels, nr_runs) << std::endl;
                  else
                    std::cout << a->median_time(
                      txt_prt, text_size, levels, nr_runs) << std::endl;
                }
              }
            }
          }
        }
      }
    }
    
    if (false) {
      remove(txt_path.c_str());
    }
  }
  return 0;
}

/******************************************************************************/
