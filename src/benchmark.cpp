/*******************************************************************************
 * src/benchmark.cpp
 *
 * Copyright (C) 2017-2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <tlx/cmdline_parser.hpp>
#include <vector>

#include <omp.h>

#include "benchmark/algorithm.hpp"
#include "util/alphabet_util.hpp"
#include "util/file_util.hpp"
#include "util/print.hpp"
#include "util/structure_decode.hpp"
#include "util/debug_assert.hpp"

#ifdef MALLOC_COUNT
#include "benchmark/malloc_count.h"
#endif // MALLOC_COUNT

#define GUARD_LOOP(x) if (!(x)) { continue; }

struct filter_by_property {
    bool exclude_other = false;
    bool exclude_self = false;

    inline bool should_keep(bool algo_property) {
        if (exclude_other && exclude_self) {
            std::cerr << "An pair of options excluded all algorithms, check "
                      << "the command line arguments" << std::endl;;
            exit(1);
        }

        if (algo_property && exclude_self) {
            return false;
        }
        if (!algo_property && exclude_other) {
            return false;
        }

        return true;
    }
};

struct {
  std::vector<std::string> file_paths;
  std::string filter_name = "";
  unsigned int word_width = 1;
  unsigned int nr_runs = 5;
  uint64_t prefix_size = 0;
  bool list_algorithms_only = false;

  filter_by_property parallel_filter;
  filter_by_property huffman_filter;
  filter_by_property matrix_filter;

  //filter_by_property external_filter;
  bool external = false;

  //filter_by_property external_in_filter;
  bool external_in = false;

  //filter_by_property external_out_filter;
  bool external_out = false;
  std::string em_dir1;
  std::string em_dir2;
  std::string em_dir3;
  std::string em_dir4;
  std::string em_dir5;
  std::string em_dir6;

  bool memory = false;
  bool check = false;
  bool debug_print = false;

  int32_t number_threads;

} global_settings;

struct {
  template <bool ext_in, bool ext_out, unsigned int width>
  int32_t run() {
    // determine input and output type of the requested algorithms
    using in_type = typename input_type<ext_in, width>::type;
    using in_type_internal = typename input_type<false, width>::type;
    using out_type = typename output_type<ext_out>::type;
    using out_type_internal = typename output_type<false>::type;
    using uint_t = typename type_for_bytes<width>::type;

    int returncode = 0;

    // get all algorithms for given types
    auto& algo_list = algorithm_list<in_type, out_type>::get_algorithm_list();
    if (global_settings.list_algorithms_only) {
      for (const auto& a : algo_list) {
        a->print_info();
      }
      return 0;
    }

    for (const auto& path : global_settings.file_paths) {
      in_type input_for_algo;
      uint64_t text_size = 0;
      uint64_t max_char = 0;
      uint64_t levels = 0;
      std::vector<uint_t> text_uint;
#ifdef MALLOC_COUNT
      malloc_count_reset_peak();
#endif
      if constexpr(!ext_in) {
        text_uint = file_to_vector<width>(path, global_settings.prefix_size);
        text_size = text_uint.size();
        max_char = reduce_alphabet(text_uint);
        levels = levels_for_max_char(max_char);
        input_for_algo = text_uint.data();
      } else {
        stxxl::syscall_file stxxl_file(path, stxxl::file::open_mode::RDONLY);
        const stxxlvector<uint_t> unreduced_vector(&stxxl_file);
        input_for_algo = stxxlvector<uint_t>();
        max_char = reduce_alphabet<uint_t>(
            unreduced_vector,
            input_for_algo,
            global_settings.prefix_size);
        text_size = input_for_algo.size();
        levels = levels_for_max_char(max_char);
      }
      for (const auto& a : algo_list) {
        GUARD_LOOP(global_settings.filter_name == "" ||
                   (a->name().compare(global_settings.filter_name) == 0));
        GUARD_LOOP(global_settings.parallel_filter.should_keep(a->is_parallel()));
        GUARD_LOOP(global_settings.huffman_filter.should_keep(a->is_huffman_shaped()));
        GUARD_LOOP(global_settings.matrix_filter.should_keep(!a->is_tree()));

        std::cout << "RESULT " << "algo=" << a->name() << ' ';
        if (global_settings.memory) {
#ifdef MALLOC_COUNT
          malloc_count_reset_peak();
          a->memory_peak(input_for_algo, text_size, levels);
          std::cout << "memory=" << malloc_count_peak() << ' ';
#else
          std::cout << "memory=no ";
#endif // MALLOC_COUNT
        }
        std::cout << "runs=" << global_settings.nr_runs << " ";
        if (global_settings.nr_runs > 0) {
          auto median_time =  a->median_time(input_for_algo, text_size, levels,
                                             global_settings.nr_runs);
          std::cout << "median_time=" << median_time << ' ';
        }
        std::cout << "input=" << path << ' '
                  << "characters=" << text_size << ' '
                  << "sigma=" << max_char + 1 << ' '
                  << "word_width=" << global_settings.word_width << ' '
                  << "threads=" << (a->is_parallel() ? global_settings.number_threads : 1)
                  << std::endl;


        if (global_settings.debug_print || global_settings.check) {
          auto ws_getter = [&](){
            if constexpr (ext_out)
              // TODO: this currently only checks structure and decoding for
              // external results (histograms should be checked, too)
              return a->compute_bitvector(input_for_algo, text_size, levels)
                  .getInternalStructure();
            else
              return a->compute_bitvector(input_for_algo, text_size, levels);
          };
          auto in_getter = [&]() {
              if constexpr (ext_in) {
                text_uint = file_to_vector<width>(
                    path, global_settings.prefix_size);
                text_size = text_uint.size();
                max_char = reduce_alphabet(text_uint);
                levels = levels_for_max_char(max_char);
                return text_uint.data();
              }
              else
                return input_for_algo;
          };
          auto structure = ws_getter();
          auto input_for_algo_check = in_getter();
          if (global_settings.debug_print) {
            print_structure(std::cout, structure, true);
          }
          if (global_settings.check) {
            if (global_settings.word_width != 1) {
              std::cout << "WARNING:"
                        << " Can only check texts over 1-byte alphabets\n";
            } else {
              auto& algo_list_check =
                  algorithm_list<in_type_internal, out_type_internal>::
                      get_algorithm_list();
              construction_algorithm<in_type_internal, out_type_internal> const* naive =
                  nullptr;

              if ((a->is_tree()) && !(a->is_huffman_shaped())) {
                naive = algo_list_check.filtered([](auto e) {
                    return e->name() == "wt_naive";
                }).at(0);
              }
              if (!(a->is_tree()) && !(a->is_huffman_shaped())) {
                naive = algo_list_check.filtered([](auto e) {
                    return e->name() == "wm_naive";
                }).at(0);
              }
              if ((a->is_tree()) && (a->is_huffman_shaped())) {
                naive = algo_list_check.filtered([](auto e) {
                    return e->name() == "wt_huff_naive";
                }).at(0);
              }
              if (!(a->is_tree()) && (a->is_huffman_shaped())) {
                naive = algo_list_check.filtered([](auto e) {
                    return e->name() == "wm_huff_naive";
                }).at(0);
              }
              CHECK(naive != nullptr);
              auto naive_wx = naive->compute_bitvector(
                  input_for_algo_check, text_size, levels);
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
                  std::cout << "Compare zeros arrays..." << std::endl;
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
                  std::cout << "Compare level " << l << "..." << std::endl;
                  auto sbs = sbvs.level_bit_size(l);
                  auto nbs = nbvs.level_bit_size(l);
                  if(check_err(sbs == nbs,
                               std::string("bit size differs on level ")
                               + std::to_string(l))) {
                    for (uint64_t bi = 0; bi < sbs; bi++) {
                      constexpr uint64_t entry_bits
                        = (sizeof(decltype(sbvs[0][0])) * 8);

                      auto sbvs_l_bi = sbvs[l][bi / entry_bits];
                      auto nbvs_l_bi = nbvs[l][bi / entry_bits];

                      if (((bi % entry_bits) == 0) && (sbvs_l_bi == nbvs_l_bi)) {
                        bi += (entry_bits - 1);
                        continue;
                      }

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
                if (!global_settings.debug_print) {
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
              if (std::equal(text_uint.begin(), text_uint.end(),
                             decoded.begin(), decoded.end())) {
                std::cout << "Output decoded OK" << std::endl;
              } else {
                std::cout << "ERROR:"
                          << "Decoded output not equal to input!"
                          << std::endl;
                std::cout << "Input:" << std::endl;
                pvec(text_uint);
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

  template <bool ext_in, bool ext_out>
  int32_t startStep2() {

    if constexpr (ext_in || ext_out) {
      std::vector<std::string> em_dirs;
      if(!global_settings.em_dir1.empty())
        em_dirs.push_back(global_settings.em_dir1);
      if(!global_settings.em_dir2.empty())
        em_dirs.push_back(global_settings.em_dir2);
      if(!global_settings.em_dir3.empty())
        em_dirs.push_back(global_settings.em_dir3);
      if(!global_settings.em_dir4.empty())
        em_dirs.push_back(global_settings.em_dir4);
      if(!global_settings.em_dir5.empty())
        em_dirs.push_back(global_settings.em_dir5);
      if(!global_settings.em_dir6.empty())
        em_dirs.push_back(global_settings.em_dir6);

      std::cout << "Setting up external memory...";
      if(em_dirs.size() > 0) {
        std::cout << std::endl;
        for(unsigned i = 0; i < em_dirs.size(); ++i) {
          std::cout << "EM buffer directory " << (i + 1) << ": " << em_dirs[i] << std::endl;
          stxxl_files::addDirectory(em_dirs[i]);
        }
      } else {
        std::cout << " No EM directories files given." << std::endl;
      }
    }

    if(global_settings.word_width == 1)
      return run<ext_in, ext_out, 1>();
    else if(global_settings.word_width == 2)
      return run<ext_in, ext_out, 2>();
    else if(global_settings.word_width == 4)
      return run<ext_in, ext_out, 4>();
    else if(global_settings.word_width == 8)
      return run<ext_in, ext_out, 8>();
    else {
      std::cerr << "You entered an invalid number of bytes per character "
                   "(parameter 'b')." << std::endl;
      return -1;
    }
  }

  template <bool ext_in>
  int32_t startStep1() {
    if(global_settings.external_out)
      return startStep2<ext_in, true>();
    else
      return startStep2<ext_in ,false>();
  }

  int32_t start() {
    if(global_settings.external_in)
      return startStep1<true>();
    else
      return startStep1<false>();
  }
} execution;

int32_t main(int32_t argc, char const* argv[]) {
  tlx::CmdlineParser cp;

  cp.set_description("Parallel Wavelet Tree and Wavelet Matrix Construction");
  cp.set_author("Florian Kurpicz <florian.kurpicz@tu-dortmund.de>\n"
                "        Marvin Löbel <loebel.marvin@gmail.com>\n"
                "        Jonas Ellert <jonas.ellert@tu-dortmund.de>\n");

  cp.add_stringlist('f', "file", global_settings.file_paths,
                    "Path(s) to the text file(s).");
  cp.add_string('n', "name", global_settings.filter_name,
                "Runs all algorithms that contain the <name> in their name");
  cp.add_uint('b', "byte", global_settings.word_width,
              "Bytes per char in the input text.");
  cp.add_uint('r', "runs", global_settings.nr_runs,
              "Number of repetitions of the construction algorithm.");
  cp.add_bytes('l', "length", global_settings.prefix_size,
               "Length of the prefix of the text that should be considered");
  cp.add_flag('\0', "list", global_settings.list_algorithms_only,
              "Print the name and description of all registered algorithms");

  cp.add_flag('\0', "parallel", global_settings.parallel_filter.exclude_other,
              "Run only parallel construction algorithms.");
  cp.add_flag('\0', "sequential", global_settings.parallel_filter.exclude_self,
              "Run only sequential construction algorithms.");

  cp.add_flag('\0', "huffman", global_settings.huffman_filter.exclude_other,
              "Run only huffman-shaped construction algorithms.");
  cp.add_flag('\0', "no_huffman", global_settings.huffman_filter.exclude_self,
              "Run only uncompressed (non-Huffman) construction algorithms");
  cp.add_flag('\0', "no_trees", global_settings.matrix_filter.exclude_other,
              "Skip all wavelet trees construction algorithms.");
  cp.add_flag('\0', "no_matrices", global_settings.matrix_filter.exclude_self,
              "Skip all wavelet matrices construction algorithms.");

  cp.add_flag('\0', "external", global_settings.external,
              "Run only external memory algorithms.");
  cp.add_flag('\0', "external_in", global_settings.external_in,
              "Run only semi-external algorithms (stream input from disk).");
  cp.add_flag('\0', "external_out", global_settings.external_out,
              "Run only semi-external algorithms (stream output to disk).");
  cp.add_string('\0', "em_dir1", global_settings.em_dir1,
                "Use specified file as external memory");
  cp.add_string('\0', "em_dir2", global_settings.em_dir2,
                "Use specified file as external memory");
  cp.add_string('\0', "em_dir3", global_settings.em_dir3,
                "Use specified file as external memory");
  cp.add_string('\0', "em_dir4", global_settings.em_dir4,
                "Use specified file as external memory");
  cp.add_string('\0', "em_dir5", global_settings.em_dir5,
                "Use specified file as external memory");
  cp.add_string('\0', "em_dir6", global_settings.em_dir6,
                "Use specified file as external memory");

  cp.add_flag('m', "memory", global_settings.memory,
              "Compute peak memory during construction.");
  cp.add_flag('c', "check", global_settings.check,
              "Check the constructed wavelet structure for validity.");
  cp.add_flag('d', "debug_print", global_settings.debug_print,
              "Output the bit vectors in a human readable format to stdout.");

  if (!cp.process(argc, argv)) {
    return -1;
  }

  global_settings.number_threads = omp_get_max_threads();

  if (global_settings.external) {
    global_settings.external_in = true;
    global_settings.external_out = true;
  }

  return execution.start();
}

/******************************************************************************/
