/*******************************************************************************
 * src/benchmark.cpp
 *
 * Copyright (C) 2017-2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <tlx/cmdline_parser.hpp>
#include <vector>

#include "benchmark/algorithm.hpp"
#include "util/alphabet_util.hpp"
#include "util/file_util.hpp"
#include "util/print.hpp"
#include "util/structure_decode.hpp"

#ifdef MALLOC_COUNT
#include "benchmark/malloc_count.h"
#endif // MALLOC_COUNT

#define GUARD_LOOP(x) if (!(x)) { continue; }

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
  std::vector<std::string> file_paths;
  std::string filter = "";
  unsigned int word_width = 1;
  unsigned int nr_runs = 5;
  uint64_t prefix_size = 0;
  bool list_algorithms_only = false;
  bool run_only_parallel = false;
  bool run_only_sequential = false;
  bool run_only_huffman = false;
  bool run_no_huffman = false;
  bool no_trees = false;
  bool no_matrices = false;
  bool memory = false;
  bool check = false;
  bool debug_print = false;

  tlx::CmdlineParser cp;

  cp.set_description("Parallel Wavelet Tree and Wavelet Matrix Construction");
  cp.set_author("Florian Kurpicz <florian.kurpicz@tu-dortmund.de>\n"
                "        Marvin Löbel <loebel.marvin@gmail.com>");

  cp.add_stringlist('f', "file", file_paths,
                    "Path(s) to the text file(s).");
  cp.add_string('n', "name", filter,
                "Runs all algorithms that contain the <name> in their name");
  cp.add_uint('b', "byte", word_width,
              "Bytes per char in the input text.");
  cp.add_uint('r', "runs", nr_runs,
              "Number of repetitions of the construction algorithm.");
  cp.add_bytes('l', "length", prefix_size,
               "Length of the prefix of the text that should be considered");
  cp.add_flag('\0', "list", list_algorithms_only,
              "Print the name and description of all registered algorithms");
  cp.add_flag('p', "parallel", run_only_parallel,
              "Run only parallel construction algorithms.");
  cp.add_flag('s', "sequential", run_only_sequential,
              "Run only sequential construction algorithms.");
  cp.add_flag('h', "huffman", run_only_huffman,
              "Run only huffman-shaped construction algorithms.");
  cp.add_flag('u', "no_huffman", run_no_huffman,
              "Run only uncompressed (non-Huffman) construction algorithms");
  cp.add_flag('\0', "no_trees", no_trees,
              "Skip all wavelet trees construction algorithms.");
  cp.add_flag('\0', "no_matrices", no_matrices,
              "Skip all wavelet matrices construction algorithms.");
  cp.add_flag('\0', "memory", memory,
              "Compute peak memory during construction.");
  cp.add_flag('c', "check", check,
              "Check the constructed wavelet structure for validity.");
  cp.add_flag('d', "debug_print", debug_print,
              "Output the bit vectors in a human readable format to stdout.");

  if (!cp.process(argc, argv)) {
    return -1;
  }

  int returncode = 0;

  auto& algo_list = algorithm_list::get_algorithm_list();
  if (list_algorithms_only) {
    for (const auto& a : algo_list) {
      a->print_info();
    }
    return 0;
  }

  for (const auto& path : file_paths) {
    std::cout << std::endl << "Text: " << path << std::endl;
    void* txt_prt = nullptr;
    uint64_t text_size = 0;
    uint64_t max_char = 0;
    uint64_t levels = 0;
    std::vector<uint8_t> text_uint8;
    std::vector<uint16_t> text_uint16;
    std::vector<uint32_t> text_uint32;
    std::vector<uint64_t> text_uint64;
#ifdef MALLOC_COUNT
    malloc_count_reset_peak();
#endif
    if (word_width == 1) {
      text_uint8 = file_to_vector<1>(path, prefix_size);
      text_size = text_uint8.size();
      max_char = reduce_alphabet(text_uint8);
      levels = levels_for_max_char(max_char);
      txt_prt = &text_uint8;
    } else if (word_width == 2) {
      text_uint16 = file_to_vector<2>(path, prefix_size);
      text_size = text_uint16.size();
      max_char = reduce_alphabet(text_uint16);
      levels = levels_for_max_char(max_char);
      txt_prt = &text_uint16;
    } else if (word_width == 4) {
      text_uint32 = file_to_vector<4>(path, prefix_size);
      text_size = text_uint32.size();
      max_char = reduce_alphabet(text_uint32);
      levels = levels_for_max_char(max_char);
      txt_prt = &text_uint32;
    } else if (word_width == 8) {
      text_uint64 = file_to_vector<8>(path, prefix_size);
      text_size = text_uint64.size();
      max_char = reduce_alphabet(text_uint64);
      levels = levels_for_max_char(max_char);
      txt_prt = &text_uint64;
    } else {
      std::cerr << "You entered an invalid number of bytes per character "
                   "(parameter 'b')." << std::endl;
      return -1;
    }
    std::cout << "Characters: " << text_size << std::endl;
#ifdef MALLOC_COUNT
    if (memory) {
      std::cout << "Memory peak text: " << malloc_count_peak() << " B, "
                << malloc_count_peak() / (1024 * 1024) << " MiB" << std::endl;
    }
#endif // MALLOC_COUNT
    for (const auto& a : algo_list) {
      GUARD_LOOP(filter == "" || (a->name().compare(filter) == 0));
      GUARD_LOOP(a->word_width() == word_width);
      GUARD_LOOP(filter_parallel(run_only_parallel, a->is_parallel()));
      GUARD_LOOP(filter_sequential(run_only_sequential, a->is_parallel()));
      GUARD_LOOP(filter_wavelet_type(a->is_tree(), no_trees, no_matrices));
      GUARD_LOOP((!run_only_huffman) || a->is_huffman_shaped());
      GUARD_LOOP((!run_no_huffman) || (!a->is_huffman_shaped()));

      std::cout << "RESULT " << "algo=" << a->name() << ' ';
      if (memory) {
#ifdef MALLOC_COUNT
        malloc_count_reset_peak();
        a->memory_peak(txt_prt, text_size, levels);
        std::cout << "memory=" << malloc_count_peak() << ' ';
#else
        std::cout << "memory=no ";
#endif // MALLOC_COUNT
      }
      std::cout << "runs=" << nr_runs << " ";
      if (nr_runs > 0) {
        std::cout << "median_time=" << a->median_time(
          txt_prt, text_size, levels, nr_runs) << ' ';
      }
      std::cout << "input=" << path << ' '
                << "characters=" << text_size << ' '
                << "sigma=" << max_char + 1 << ' '
                << "word_width=" << word_width << std::endl;

      if (debug_print || check) {
        auto structure =
          a->compute_bitvector(txt_prt, text_size, levels);
        if (debug_print) {
          print_structure(std::cout, structure, true);
        }
        if (check) {
          if (word_width != 1) {
            std::cout << "WARNING:"
                      << " Can only check texts over 1-byte alphabets\n";
          } else {
            construction_algorithm const* naive = nullptr;
            if ((a->is_tree()) && !(a->is_huffman_shaped())) {
              naive = algo_list.filtered([](auto e) {
                return e->name() == "wt_naive";
              }).at(0);
            }
            if (!(a->is_tree()) && !(a->is_huffman_shaped())) {
              naive = algo_list.filtered([](auto e) {
                return e->name() == "wm_naive";
              }).at(0);
            }
            if ((a->is_tree()) && (a->is_huffman_shaped())) {
              naive = algo_list.filtered([](auto e) {
                return e->name() == "wt_huff_naive";
              }).at(0);
            }
            if (!(a->is_tree()) && (a->is_huffman_shaped())) {
              naive = algo_list.filtered([](auto e) {
                return e->name() == "wm_huff_naive";
              }).at(0);
            }
            assert(naive != nullptr);
            auto naive_wx =
            naive->compute_bitvector(txt_prt, text_size, levels);
            bool err_trigger = false;
            auto check_err = [&](bool cond, auto const& msg) {
              if (!cond) {
                std::cout << "ERROR: " << msg << std::endl;
                err_trigger = true;
              }
              return cond;
            };

            if (check_err(structure.levels() == naive_wx.levels(),
                          "structures have different level sizes")) {
              if (!a->is_tree()) {
                size_t sl = structure.levels();
                check_err(structure.zeros().size() == sl,
                          "structure zeros too short");
                if (sl > 0) {
                  auto sz = structure.zeros();
                  auto nz = naive_wx.zeros();
                  sz.pop_back();
                  nz.pop_back();
                  check_err(sz == nz, "zeros arrays differ");
                }
              }
              auto& sbvs = structure.bvs();
              auto& nbvs = naive_wx.bvs();
              for (size_t l = 0; l < structure.levels(); l++) {
                auto sbs = sbvs.level_bit_size(l);
                auto nbs = nbvs.level_bit_size(l);
                if(check_err(sbs == nbs,
                             std::string("bit size differs on level ")
                             + std::to_string(l))) {
                  for (uint64_t bi = 0; bi < sbs; bi++) {
                    if(!check_err(
                      bit_at(sbvs[l], bi) == bit_at(nbvs[l], bi),
                       std::string("bit ")
                       + std::to_string(bi)
                       + " differs on level "
                       + std::to_string(l))) {
                      break;
                    }
                  }
                }
              }
            }
            if (err_trigger) {
              returncode = -2;
            } else {
              std::cout << "Output structurally OK" << std::endl;
            }

            if (err_trigger) {
              if (!debug_print) {
                std::cout << "Output:\n";
                print_structure(std::cout, structure, true);
              }
              std::cout << "Naive result as comparison:\n";
              print_structure(std::cout, naive_wx, true);
            }

            auto pvec = [](auto const& v) {
              std::cout << "[";
              for (auto e : v) {
                  std::cout << uint64_t(e) << ", ";
              }
              std::cout << "]\n";
            };

            std::string decoded_ = decode_structure(structure);
            std::vector<uint8_t> decoded(decoded_.begin(), decoded_.end());
            if (std::equal(text_uint8.begin(), text_uint8.end(),
                           decoded.begin(), decoded.end())) {
              std::cout << "Output decoded OK" << std::endl;
            } else {
              std::cout << "ERROR:"
                        << "Decoded output not equal to input!"
                        << std::endl;
              std::cout << "Input:" << std::endl;
              pvec(text_uint8);
              std::cout << "Decoded:" << std::endl;
              pvec(decoded);
            }
          }
          std::cout << std::endl;
        }
      }
    }
  }
  return returncode;
}

/******************************************************************************/
