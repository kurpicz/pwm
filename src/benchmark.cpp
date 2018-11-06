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

#include "benchmark/algorithm.hpp"
#include "util/alphabet_util.hpp"
#include "util/file_util.hpp"
#include "util/print.hpp"
#include "util/structure_decode.hpp"

#ifdef MALLOC_COUNT
#include "benchmark/malloc_count.h"
#endif // MALLOC_COUNT

#define GUARD_LOOP(x) if (!(x)) { continue; }

struct {
  std::vector<std::string> file_paths;
  std::string filter = "";
  unsigned int word_width = 1;
  unsigned int nr_runs = 5;
  uint64_t prefix_size = 0;
  bool list_algorithms_only = false;
  bool run_only_parallel = false;
  bool run_only_sequential = false;
  bool run_only_huffman = false;
  bool run_only_no_huffman = false;
  bool no_trees = false;
  bool no_matrices = false;

  bool external = false;
  bool external_in = false;
  bool external_out = false;
  std::string em_file1;
  std::string em_file2;
  std::string em_file3;
  std::string em_file4;
  std::string em_file5;
  std::string em_file6;

  bool memory = false;
  bool check = false;
  bool debug_print = false;

  auto filter_parallel(bool is_parallel) {
    return (!run_only_parallel || is_parallel);
  }

  auto filter_sequential(bool is_parallel) {
    return (!run_only_sequential || !is_parallel);
  }

  auto filter_wavelet_type(bool is_tree) {
    return (is_tree ? !no_trees : !no_matrices);
  }

  auto filter_huffman(bool is_huffman) {
    return (!run_only_huffman || is_huffman);
  }

  auto filter_no_huffman(bool is_huffman) {
    return (!run_only_no_huffman || !is_huffman);
  }

} global_settings;

struct {
  template <bool ext_in, bool ext_out, unsigned int width>
  int32_t run() {
    // determine input and output type of the requested algorithms
    using in_type = typename input_type<ext_in, width>::type;
    using out_type = typename output_type<ext_out>::type;
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
      std::cout << std::endl << "Text: " << path << std::endl;
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
        text_size = std::min(uint64_t(unreduced_vector.size()), global_settings.prefix_size);
        input_for_algo = stxxlvector<uint_t>();
        max_char = reduce_alphabet<uint_t>(unreduced_vector, input_for_algo);
        levels = levels_for_max_char(max_char);
      }
      std::cout << "Characters: " << text_size << std::endl;
#ifdef MALLOC_COUNT
      if (global_settings.memory) {
      std::cout << "Memory peak text: " << malloc_count_peak() << " B, "
                << malloc_count_peak() / (1024 * 1024) << " MiB" << std::endl;
    }
#endif // MALLOC_COUNT
      for (const auto& a : algo_list) {
        GUARD_LOOP(global_settings.filter == "" || (a->name().compare(global_settings.filter) == 0));
        GUARD_LOOP(global_settings.filter_parallel(a->is_parallel()));
        GUARD_LOOP(global_settings.filter_sequential(a->is_parallel()));
        GUARD_LOOP(global_settings.filter_wavelet_type(a->is_tree()));
        GUARD_LOOP(global_settings.filter_huffman(a->is_huffman_shaped()));
        GUARD_LOOP(global_settings.filter_no_huffman(a->is_huffman_shaped()));

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
          std::cout << "median_time=" << a->median_time(
              input_for_algo, text_size, levels, global_settings.nr_runs) << ' ';
        }
        std::cout << "input=" << path << ' '
                  << "characters=" << text_size << ' '
                  << "sigma=" << max_char + 1 << ' '
                  << "word_width=" << global_settings.word_width << std::endl;

        if constexpr (!ext_in && !ext_out){
          if (global_settings.debug_print || global_settings.check) {
            auto structure =
                a->compute_bitvector(input_for_algo, text_size, levels);
            if (global_settings.debug_print) {
              print_structure(std::cout, structure, true);
            }
            if (global_settings.check) {
              if (global_settings.word_width != 1) {
                std::cout << "WARNING:"
                          << " Can only check texts over 1-byte alphabets\n";
              } else {
                construction_algorithm<in_type, out_type> const* naive = nullptr;
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
                    naive->compute_bitvector(input_for_algo, text_size, levels);
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
    }
    return returncode;
  }

  template <bool ext_in, bool ext_out>
  int32_t startStep2() {

    if constexpr (ext_in || ext_out) {
      std::vector<std::string> em_files;
      if(!global_settings.em_file1.empty())
        em_files.push_back(global_settings.em_file1);
      if(!global_settings.em_file2.empty())
        em_files.push_back(global_settings.em_file2);
      if(!global_settings.em_file3.empty())
        em_files.push_back(global_settings.em_file3);
      if(!global_settings.em_file4.empty())
        em_files.push_back(global_settings.em_file4);
      if(!global_settings.em_file5.empty())
        em_files.push_back(global_settings.em_file5);
      if(!global_settings.em_file6.empty())
        em_files.push_back(global_settings.em_file6);

      std::cout << "Setting up external memory...";
      if(em_files.size() > 0) {
        std::cout << std::endl;
        for(unsigned i = 0; i < em_files.size(); ++i) {
          std::cout << "EM buffer file " << (i + 1) << ": " << em_files[i] << std::endl;
          stxxl_files::addFile(em_files[i]);
        }
      } else {
        std::cout << " No EM buffer files given." << std::endl;
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
  cp.add_string('n', "name", global_settings.filter,
                "Runs all algorithms that contain the <name> in their name");
  cp.add_uint('b', "byte", global_settings.word_width,
              "Bytes per char in the input text.");
  cp.add_uint('r', "runs", global_settings.nr_runs,
              "Number of repetitions of the construction algorithm.");
  cp.add_bytes('l', "length", global_settings.prefix_size,
               "Length of the prefix of the text that should be considered");
  cp.add_flag('\0', "list", global_settings.list_algorithms_only,
              "Print the name and description of all registered algorithms");
  cp.add_flag('p', "parallel", global_settings.run_only_parallel,
              "Run only parallel construction algorithms.");
  cp.add_flag('s', "sequential", global_settings.run_only_sequential,
              "Run only sequential construction algorithms.");
  cp.add_flag('h', "huffman", global_settings.run_only_huffman,
              "Run only huffman-shaped construction algorithms.");
  cp.add_flag('u', "no_huffman", global_settings.run_only_no_huffman,
              "Run only uncompressed (non-Huffman) construction algorithms");
  cp.add_flag('\0', "no_trees", global_settings.no_trees,
              "Skip all wavelet trees construction algorithms.");
  cp.add_flag('\0', "no_matrices", global_settings.no_matrices,
              "Skip all wavelet matrices construction algorithms.");

  cp.add_flag('\0', "external", global_settings.external,
              "Run only external memory algorithms.");
  cp.add_flag('\0', "external_in", global_settings.external_in,
              "Run only semi-external algorithms (stream input from disk).");
  cp.add_flag('\0', "external_out", global_settings.external_out,
              "Run only semi-external algorithms (stream output to disk).");
  cp.add_string('\0', "em_file1", global_settings.em_file1,
                "Use specified file as external memory");
  cp.add_string('\0', "em_file2", global_settings.em_file2,
                "Use specified file as external memory");
  cp.add_string('\0', "em_file3", global_settings.em_file3,
                "Use specified file as external memory");
  cp.add_string('\0', "em_file4", global_settings.em_file4,
                "Use specified file as external memory");
  cp.add_string('\0', "em_file5", global_settings.em_file5,
                "Use specified file as external memory");
  cp.add_string('\0', "em_file6", global_settings.em_file6,
                "Use specified file as external memory");

  cp.add_flag('\0', "memory", global_settings.memory,
              "Compute peak memory during construction.");
  cp.add_flag('c', "check", global_settings.check,
              "Check the constructed wavelet structure for validity.");
  cp.add_flag('d', "debug_print", global_settings.debug_print,
              "Output the bit vectors in a human readable format to stdout.");

  if (!cp.process(argc, argv)) {
    return -1;
  }

  if (global_settings.external) {
    global_settings.external_in = true;
    global_settings.external_out = true;
  }

  return execution.start();
}

/******************************************************************************/
