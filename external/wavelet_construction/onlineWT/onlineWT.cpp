/*******************************************************************************
 * onlineWT/onlineWT.cpp
 *
 * Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <experimental/filesystem> // compile with gcc 7.3 because of cilk
#include <iostream>
#include <vector>

#include "src/bytearray.hpp"
#include "src/cstringutil.hpp"
#include "src/strstream.hpp"
#include "src/wavtree.hpp"

#include <tlx/cmdline_parser.hpp>

#ifdef MALLOC_COUNT
#include "benchmark/malloc_count.h"
#elif defined GET_RUSAGE
#include "../getmemory.hpp"
#endif // MALLOC_COUNT

int32_t main(int32_t argc, char* argv[]) {
  std::string iFile;
  int rounds = 5;
  uint64_t prefix_size = 0;

  tlx::CmdlineParser cp;
  cp.add_param_string("input", iFile, "Path to the input text");
  cp.add_int('r', "rounds", rounds, "Number of executions of the algorithm");
  cp.add_bytes('l', "length", prefix_size,
                "Length of the prefix of the text that should be considered");

  if (!cp.process(argc, argv)) {
    return -1;
  }

  namespace fs = std::experimental::filesystem;

  fs::path in_path = iFile.c_str();
  fs::path tmp_path = iFile + "_wx_experiments_tmp";
  std::string tmp_path_str = tmp_path.string();

  fs::copy_file(in_path, tmp_path);

  if (prefix_size > 0 && prefix_size < fs::file_size(tmp_path)) {
    fs::resize_file(tmp_path, prefix_size);
  }

  // Here, we only compute the alphabet for our statistics (it is not used
  // during the computation of the wavelet tree.
  strstream *sst = strstream_open_file(tmp_path_str.data());
  byte_t *ab_table = bytearr_new(256);
  size_t n=0, m=0;
  for (int c=0; (c=(strstream_getc(sst)))!=EOF; n++) {
    if (!ab_table[c]) {
      ab_table[c]=1; 
      m++;
    }
  }
  char *ab_str = cstr_new(m);
  size_t k = 0;
  for (int i=0; i<256; i++)
    if (ab_table[i]==1) ab_str[k++] = (char)i;

  alphabet *ab = alphabet_new(m, ab_str);

  std::cout << "RESULT algo=wt_online ";

  wavtree *wt = nullptr;

#ifdef MALLOC_COUNT
  malloc_count_reset_peak();
  wt = wavtree_new_online_from_stream(sst, WT_BALANCED);
  std::cout << "memory=" << malloc_count_peak() << ' ';
#elif defined GET_RUSAGE
  wt = wavtree_new_online_from_stream(sst, WT_BALANCED);
  std::cout << "memory=" << getPeakRSS() << ' ';
#else
  std::cout << "memory=no ";
#endif
  
  std::cout << "runs=" << rounds << ' ';
  std::vector<float> times;
  for (int i=0; i < rounds; i++) {
    free(wt);
    auto begin_time = std::chrono::high_resolution_clock::now();
    wt = wavtree_new_online_from_stream(sst, WT_BALANCED);
    auto end_time = std::chrono::high_resolution_clock::now();
    times.emplace_back(static_cast<float>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - begin_time).count()));
  }
  std::sort(times.begin(), times.end());
  std::cout << "median_time=";
  if (rounds % 2 == 0) {
    std::cout << (times[rounds >> 1] + times[(rounds >> 1) - 1]) / 2;
  } else {
    std::cout << times[rounds >> 1]; 
  }
  std::cout << " input=" << iFile << ' '
            << "characters=" << n << ' '
            << "sigma=" << ab_size(ab) << ' '
            << "word_width=1" << ' '
            << "threads=1" << std::endl;

  fs::remove(tmp_path);

  return 0;
}

/******************************************************************************/
