/*******************************************************************************
 * benchmark/construct_wt.cpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <wt_pc.hpp>
#include <wt_ppc.hpp>
#include <wt_ps.hpp>
#include <wt_pps.hpp>

template <typename AlphabetType>
void ConstructWT(std::vector<AlphabetType>& text, const bool already_reduced) {

  uint64_t max_char = 0;
  if (!already_reduced) {
    std::unordered_map<AlphabetType, uint64_t> word_list;
    for (const AlphabetType& character : text) {
      auto result = word_list.find(character);
      if (result == word_list.end()) {
        word_list.emplace(character, max_char++);
      }
    }
    --max_char;
    for (uint64_t i = 0; i < text.size(); ++i) {
      text[i] = static_cast<AlphabetType>(word_list.find(text[i])->second);
    }
  } else {
    for (const AlphabetType& character : text) {
      if (character > max_char) {
        max_char = character;
      }
    }
  }
  std::cout << "Greatest character: " << static_cast<uint64_t>(max_char)
            << std::endl;

  AlphabetType levels = 0;
  while (max_char) {
    max_char >>= 1;
    ++levels;
  }
  std::cout << "Levels in WM: " << static_cast<uint64_t>(levels) << std::endl;

#ifdef TIMING // Get the average construction time over 5 (default) runs.
  auto t1 = std::chrono::high_resolution_clock::now();
  for (size_t run = 0; run < RUNS; ++run) {
    WT_TYPE<AlphabetType, uint32_t> wt(text, text.size(), levels);
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto seq_time_sorting =
    ((std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count())
      / RUNS);
  std::cout << "Construction time: "
            << static_cast<float>(seq_time_sorting) / 1000 << " seconds."
            << std::endl;
#elif CHECK // Check the correctness of the construction algorithm.
  std::vector<uint64_t*> wm_bv;
  std::vector<uint32_t> wm_zeros;
  std::vector<uint64_t*> wm_naive_bv;
  std::vector<uint32_t> wm_naive_zeros;

  WT_TYPE<AlphabetType, uint32_t> wt(text, text.size(), levels);
  // std::tie(wm_bv, wm_zeros) = wt.get_bv();

  // wt_naive<AlphabetType, uint32_t> wm_naive(text, text.size(), levels);
  // std::tie(wm_naive_bv, wm_naive_zeros) = wm_naive.get_bv_and_zeros();

  // for (AlphabetType level = 0; level < levels; ++level) {
  //   for (uint64_t i = 0; i < (text.size() >> 6); ++i) {
  //     if (wm_bv[level][i] != wm_naive_bv[level][i]) {
  //       std::cout << "Error in level " << static_cast<uint64_t>(level)
  //                 << " at position " << i << std::endl;
  //       std::exit(EXIT_FAILURE);
  //     }
  //   }
  //   if (wm_zeros[level] != wm_naive_zeros[level]) {
  //     std::cout << "Zeros in level " << static_cast<uint64_t>(level)
  //               << " not matching." << std::endl;
  //     std::cout << wm_zeros[level] << " given, while " << wm_naive_zeros[level]
  //               << " expected." << std::endl;
  //     std::exit(EXIT_FAILURE);
  //   }
  // }
  std::cout << "Algorithm working corretly." << std::endl;
#elif MEMORY // Measure the memory consumption of the construction algorithm.
  WT_TYPE<AlphabetType, uint32_t> wt(text, text.size(), levels);
#endif
}

int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cout << argv[0] << " requires at least two parameter: "
              << argv[0] << " [file name] [symbol width in byte]"
                         << " [already reduced (optional, any sequence will"
                         << " turn this option on!)]" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  uint32_t byte_width = std::stoi(std::string(argv[2]));
  if (byte_width > 8 || 
    (byte_width != 1 && ((byte_width == 6) || ((byte_width % 2) == 1)))) {
    std::cout << "[symbol with in bytes] must be 1, 2, 4 or 8." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const bool already_reduced = (argc >= 4);

  std::cout << argv[0] << " " << argv[1];
  if (already_reduced) {
    std::cout << " with already reduced alphabet." << std::endl;
  } else {
    std::cout << " with not yet reduced alphabet." << std::endl;
  }

  std::ifstream stream(argv[1], std::ios::in | std::ios::binary);
  stream.seekg(0, std::ios::end);
  uint64_t size = stream.tellg();
  stream.seekg(0);

  if (byte_width == 1) {
    std::vector<uint8_t> text(size);
    stream.read(reinterpret_cast<char*>(text.data()), size);
    ConstructWT(text, already_reduced);
  } else if (byte_width == 2) {
    std::vector<uint16_t> text(size >> 1);
    stream.read(reinterpret_cast<char*>(text.data()), size);
    ConstructWT(text, already_reduced);
  } else if (byte_width == 4) {
    std::vector<uint32_t> text(size >> 2);
    stream.read(reinterpret_cast<char*>(text.data()), size);
    ConstructWT(text, already_reduced);
  } else if (byte_width == 8) {
    std::vector<uint64_t> text(size >> 3);
    stream.read(reinterpret_cast<char*>(text.data()), size);
    ConstructWT(text, already_reduced);
  }
  stream.close();

  return 0;
}

/******************************************************************************/