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
#include "benchmark/file_util.hpp"

auto filter_parallel(bool only_parallel, bool is_parallel) {
  return (!only_parallel || is_parallel);
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
    "Number of repetitions of the construction algorithm.", false, 5, "int32_t");
  TCLAP::SwitchArg run_only_parallel_arg("p", "parallel",
    "Run only parallel construction algorithms.", false);
  cmd.add(run_only_parallel_arg);
  TCLAP::SwitchArg no_trees_arg("m", "no_trees",
    "Skip all wavelet trees construction algorithms.", false);
  cmd.add(no_trees_arg);
  TCLAP::SwitchArg no_matrices_arg("t", "no_matrices",
    "Skip all wavelet matrices construction algorithms.", false);
  cmd.add(no_matrices_arg);
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
  const bool no_trees = no_trees_arg.getValue();
  const bool no_matrices = no_matrices_arg.getValue();

  for (const auto& path : file_paths) {
    auto text = file_to_vector<1>(path);
    uint64_t levels = reduce_alphabet(text);
    for (const auto& a : algo_list) {
      if (filter == "" || (a->name().find(filter) != std::string::npos)) {
        if (a->word_width() == word_width) {
          if (filter_parallel(run_only_parallel, a->is_parallel())) {
            if (filter_wavelet_type(a->is_tree(), no_trees, no_matrices)) {
              a->print_info();
              std::cout << a->median_time(&text, text.size(), levels, nr_runs)
                        << std::endl;
            }
          }
        }
      }
    }
  }
  return 0;
}

/******************************************************************************/
