/*******************************************************************************
 * src/benchmark.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <tclap/CmdLine.h>
#include <vector>

#include "benchmark/algorithm.hpp"
#include "util/alphabet_util.hpp"
#include "util/file_util.hpp"

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
  cmd.parse( argc, argv );

  auto& algo_list = algorithm_list::get_algorithm_list();
  if (list_all_algorithms.getValue()) {
    for (const auto& a : algo_list) {
      a->print_info();
    }
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

  for (const auto& path : file_paths) {
    std::cout << std::endl << "Text: " << path << std::endl;
    void* txt_prt = nullptr;
    uint64_t text_size = 0;
    uint64_t levels = 0;
    std::vector<uint8_t> text_uint8;
    std::vector<uint16_t> text_uint16;
    std::vector<uint32_t> text_uint32;
    std::vector<uint64_t> text_uint64;
#ifdef MALLOC_COUNT
    malloc_count_reset_peak();
#endif
    if (word_width == 1) {
      text_uint8 = file_to_vector<1>(path);
      text_size = text_uint8.size();
      levels = reduce_alphabet(text_uint8);
      txt_prt = &text_uint8;
    } else if (word_width == 2) {
      text_uint16 = file_to_vector<2>(path);
      text_size = text_uint16.size();
      levels = reduce_alphabet(text_uint16);
      txt_prt = &text_uint16;
    } else if (word_width == 4) {
      text_uint32 = file_to_vector<4>(path);
      text_size = text_uint32.size();
      levels = reduce_alphabet(text_uint32);
      txt_prt = &text_uint32;
    } else if (word_width == 8) {
      text_uint64 = file_to_vector<8>(path);
      text_size = text_uint64.size();
      levels = reduce_alphabet(text_uint64);
      txt_prt = &text_uint64;
    } else {
      std::cerr << "You entered an invalid number of bytes per character "
                   "(parameter 'b')." << std::endl;
      return -1;
    }
    std::cout << "Characters: " << text_size << std::endl;
#ifdef MALLOC_COUNT
    std::cout << "Memory peak text: " << malloc_count_peak() << ", MB: "
              << malloc_count_peak() / (1024 * 1024) << std::endl;
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
  }
  return 0;
}

/******************************************************************************/
