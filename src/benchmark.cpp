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
#include <sstream>

#include <omp.h>

#include "benchmark/algorithm.hpp"
#include "util/alphabet_util.hpp"
#include "util/file_util.hpp"
#include "util/print.hpp"
#include "util/structure_decode.hpp"
#include "util/stats.hpp"
#include "util/debug_assert.hpp"
#include "util/meminfo.hpp"
#include "external_memory/internal_memory_bound.hpp"

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
  std::string filter_name_contains = "";
  std::string filter_name_not_contains = "";
  unsigned int word_width = 1;
  unsigned int nr_runs = 5;
  uint64_t prefix_size = 0;
  uint64_t memory = 0;
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
  std::string em_dirs;

  bool diskbench = false;

  bool check = false;
  bool debug_print = false;

  int32_t number_threads;

  bool name_filter (const std::string& name) const {
    return ((filter_name == "") ||
            (name.compare(filter_name) == 0)) &&
        ((filter_name_contains == "") ||
            (name.find(filter_name_contains) != std::string::npos)) &&
        ((filter_name_not_contains == "") ||
            (name.find(filter_name_not_contains) == std::string::npos));
  }

} global_settings;

struct {
  template <bool ext_in, bool ext_out, unsigned int width>
  int32_t run() {
    // determine input and output type of the requested algorithms
    using in_type = typename input_type<ext_in, width>::type;
    //using out_type = typename output_type<ext_out>::type;
    using uint_t = typename type_for_bytes<width>::type;


    int returncode = 0;

    // get all algorithms for given types
    auto& algo_list =
        algorithm_list<ext_in, ext_out, width>::get_algorithm_list();
    if (global_settings.list_algorithms_only) {
      for (const auto& a : algo_list) {
        a->print_info();
      }
      return 0;
    }

    for (const auto& path : global_settings.file_paths) {
      //std::cout << std::endl << "Text: " << path << std::endl;
      #ifdef MALLOC_COUNT
      malloc_count_reset_peak();
      uint64_t malloc_count_base = malloc_count_peak();
      #endif

      in_type input_for_algo;
      uint64_t text_size = 0;
      uint64_t max_char = 0;
      uint64_t levels = 0;
      std::vector<uint_t> text_uint;
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


      // std::cout << "Characters: " << text_size << std::endl;
      #ifdef MALLOC_COUNT
      malloc_count_reset_peak();
      uint64_t malloc_count_text = malloc_count_peak() - malloc_count_base;
      //std::cout << "Memory peak text: " << malloc_count_text << " B, "
      //          << malloc_count_text / (1024 * 1024) << " MiB" << std::endl;
      #endif // MALLOC_COUNT

      for (const auto& a : algo_list) {
        GUARD_LOOP(global_settings.name_filter(a->name()));
        GUARD_LOOP(global_settings.parallel_filter.should_keep(a->is_parallel()));
        GUARD_LOOP(global_settings.huffman_filter.should_keep(a->is_huffman_shaped()));
        GUARD_LOOP(global_settings.matrix_filter.should_keep(!a->is_tree()));

        std::cout << "RESULT " << "algo=" << a->name() << ' ';

        std::cout << "runs=" << global_settings.nr_runs << " " << std::flush;
        if (global_settings.nr_runs > 0) {
          auto stats_object = a->median_time_stats(
              input_for_algo, text_size, levels, global_settings.nr_runs);
          std::cout << stats_object << ' ';
          #ifdef MALLOC_COUNT
          std::cout << "memory_text=" << malloc_count_text << ' ';
          std::cout << "memory=" << stats_object.get_total_memory() + malloc_count_text << ' ';
          #endif // MALLOC_COUNT
        }

        auto bit_size = a->huffman_bit_size(input_for_algo, text_size, levels);

        std::cout << "input=" << path << ' '
                  << "characters=" << text_size << ' '
                  << "sigma=" << max_char + 1 << ' '
                  << "word_width=" << global_settings.word_width << ' ';
        if(a->is_huffman_shaped()) {
          std::cout << "bit_size=" << bit_size << ' ';
        }
        std::cout << "threads=" << (a->is_parallel() ? global_settings.number_threads : 1)
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
                  algorithm_list<false, false, width>::
                      get_algorithm_list();
              construction_algorithm<false, false, width> const* naive =
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
      std::cout << "Setting up external memory...";
      if(global_settings.em_dirs.size() > 0) {
        std::cout << std::endl;
        std::stringstream em_dirs_stream(global_settings.em_dirs);
        std::string em_dir;
        unsigned i = 0;
        while(std::getline(em_dirs_stream, em_dir, ':')) {
          std::cout << "EM buffer directory " << (i += 1) << ": " << em_dir << std::endl;
          stxxl_files::addDirectory(em_dir);
        }
      } else {
        std::cout << " No EM directories files given." << std::endl;
        volatile stxxlvector<bool> allocate;
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
                "Runs the named algorithm.");
  cp.add_string('\0', "contains", global_settings.filter_name_contains,
                "Runs all algorithms that contain the <name> in their name");
  cp.add_string('\0', "not_contains", global_settings.filter_name_not_contains,
                "Runs all algorithms that do not contain the <name> in their name");
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
  cp.add_bytes('\0', "memory", global_settings.memory,
              "Maximum internal memory for external memory algorithms.");

  cp.add_flag('\0', "external_in", global_settings.external_in,
              "Run only semi-external algorithms (stream input from disk).");
  cp.add_flag('\0', "external_out", global_settings.external_out,
              "Run only semi-external algorithms (stream output to disk).");
  cp.add_string('\0', "em_dirs", global_settings.em_dirs,
                    "Use the given directories as external memory"
                    "(split multiple directories using ':').");

  cp.add_flag('\0', "diskbench", global_settings.diskbench,
              "Benchmark the current STXXL configuration.");

  cp.add_flag('c', "check", global_settings.check,
              "Check the constructed wavelet structure for validity.");
  cp.add_flag('d', "debug_print", global_settings.debug_print,
              "Output the bit vectors in a human readable format to stdout.");

  if (!cp.process(argc, argv)) {
    return -1;
  }

  #pragma omp parallel
  {
    #pragma omp single
    global_settings.number_threads = omp_get_max_threads();
  }

  if (global_settings.external) {
    global_settings.external_in = true;
    global_settings.external_out = true;
  }
  else if (global_settings.external_in && global_settings.external_out) {
    global_settings.external = true;
  }

  if (global_settings.external || global_settings.diskbench) {
    const uint64_t total_mem = 1024ULL * get_mem_total();
    const uint64_t avail_mem = 1024ULL * get_mem_available();

    const uint64_t total_mem_mib = total_mem / (1024ULL * 1024);
    const uint64_t avail_mem_mib = avail_mem / (1024ULL * 1024);;

    if (total_mem == 0) {
      std::cerr << "Total system memory:     "
                << " Could not read /proc/meminfo.";
    } else {
      std::cout << "Total system memory:     "
                << total_mem_mib << "MiB" << std::endl;
    }

    if (avail_mem == 0) {
      std::cerr << "Available system memory: "
                << " Could not read /proc/meminfo.";
    } else {
      std::cout << "Available system memory: "
                << avail_mem_mib << "MiB" << std::endl;
    }

    std::cout << "Available threads:       "
              << global_settings.number_threads << std::endl;
    std::cout << std::endl;
    if (global_settings.memory == 0) {
      std::cout << "No maximum internal memory usage "
                << "specified (argument --memory)." << std::endl;
      if (avail_mem > 0) {
        global_settings.memory = (avail_mem * 9) / 10;
        const uint64_t mib = global_settings.memory / (1024ULL * 1024);
        std::cout << "Using at most " << mib << "MiB of memory, "
                  << "which is 90% of the available memory." << std::endl;
      }
      else {
        global_settings.memory = 4ULL * 1024 * 1024 * 1024;
        const uint64_t mib = global_settings.memory / (1024ULL * 1024);
        std::cout << "Using at most " << mib << "MiB of memory." << std::endl;
      }
    }
    else {
      std::cout << "Memory limit: "
                << global_settings.memory / (1024ULL * 1024)
                << "MiB" << std::endl;
    }
    internal_memory_bound::value() = global_settings.memory;
    std::cout << std::endl;
  }

  if (global_settings.diskbench) {
    const uint64_t bytes = (global_settings.memory / 8) * 8;
    const uint64_t words = bytes / 8;

    const uint64_t runs =
        std::max((unsigned int)(1), global_settings.nr_runs);

    std::vector<statistics<true>> stats_w(runs);
    std::vector<statistics<true>> stats_r(runs);

    stxxlvector<uint64_t> vec_read, vec_write;
    vec_read.resize(words * runs);
    vec_write.resize(words * runs);

    std::cout << "Disk bench: " << runs << " repetitions of "
              << (bytes / 1024ULL / 1024) << "MiB." << std::endl;
    std::cout << "Filling EM vector.." << std::endl;
    {
      stxxlwriter<uint64_t> writer(vec_read.begin());
      for (uint64_t run = 0; run < runs; ++run) {
        for (uint64_t i = 0; i < words; ++i) {
          writer << i;
        }
      }
    }

    stxxlreader<uint64_t> read_bench(vec_read);
    stxxlwriter<uint64_t> write_bench(vec_write.begin());

    for (uint64_t run = 0; run < runs; ++run) {
      uint64_t * mem = (uint64_t *) malloc(bytes);

      std::cout << "Starting read test " << (run + 1) << "..." << std::endl;
      {
        stats_r[run].start();
        for (uint64_t i = 0; i < words; ++i) {
          read_bench >> mem[i];
        }
        stats_r[run].finish();
      }
      std::cout << "Done." << std::endl;

      std::cout << "Starting write test " << (run + 1) << "..." << std::endl;
      {
        stats_w[run].start();
        for (uint64_t i = 0; i < words; ++i) {
          write_bench << mem[i];
        }
        stats_w[run].finish();
      }
      std::cout << "Done." << std::endl;

      delete mem;
    }

    std::sort(stats_r.begin(), stats_r.end());
    std::sort(stats_w.begin(), stats_w.end());

    std::cout << "RESULT algo=stxxlbench_read characters=" << bytes
              << " mibs=" << ((1000.0 * bytes) / stats_r[(runs - 1) >> 1].get_total_time()) / 1024 / 1024
              << " " << stats_r[(runs - 1) >> 1] << " word_width=1" << std::endl;

    std::cout << "RESULT algo=stxxlbench_write characters=" << bytes
              << " mibs=" << ((1000.0 * bytes) / stats_w[(runs - 1) >> 1].get_total_time()) / 1024 / 1024
              << " " << stats_w[(runs - 1) >> 1] << " word_width=1" << std::endl;

    return 0;
  }
  else return execution.start();
}

/******************************************************************************/
