/*******************************************************************************
 * include/wx_huff_naive.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>

#include "arrays/bit_vectors.hpp"
#include "construction/wavelet_structure.hpp"
#include "huffman/huff_codes.hpp"

// TODO: temporary
#include "util/print.hpp"
#include <unordered_map>

template <typename AlphabetType, bool is_tree_>
class wx_huff_naive {

public:
  static constexpr bool    is_parallel = false;
  static constexpr bool    is_tree     = is_tree_;
  static constexpr uint8_t word_width  = sizeof(AlphabetType);
  static constexpr bool    is_huffman_shaped = true;

  static wavelet_structure compute(AlphabetType const* const text,
    const uint64_t size, const uint64_t /*levels*/) {

    canonical_huff_codes<AlphabetType, is_tree> codes(text, size);

    if(size == 0) {
      if constexpr (is_tree) {
        return wavelet_structure_tree_huffman<AlphabetType>(std::move(codes));
      } else {
        return wavelet_structure_matrix_huffman<AlphabetType>(std::move(codes));
      }
    }

    // [snapshot]: got text, size, created huff codes
    const uint64_t levels = codes.levels();

    auto _bv = huff_bit_vectors(levels, codes.level_sizes());
    auto& bv = _bv.raw_data();
    auto _zeros = std::vector<size_t>();
    if constexpr (!is_tree) {
        _zeros = std::vector<size_t>(levels, 0);
    }

    std::vector<AlphabetType> local_text(size);
    for(size_t i = 0; i < size; i++) {
      local_text[i] = text[i];
    }

    // Construct each level top-down
    for (uint64_t level = 0; level < levels; ++level) {
      std::vector<AlphabetType> text0;
      std::vector<AlphabetType> text1;

      // Insert the level-th MSB in the bit vector of the level (in text order)
      uint32_t cur_pos = 0;
      // Eat prefix in whole 64bit word chunks.
      // Modification of local_text ensures its only as long as the level.
      for (; cur_pos + 64 <= local_text.size(); cur_pos += 64) {
        /*
        std::cout << "level " << level
                  << " cur_pos " << cur_pos
                  << " size " << local_text.size()
                  << " ok\n";
        */
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < 64; ++i) {
          const code_pair cp = codes.encode_symbol(local_text[cur_pos + i]);
          word <<= 1;
          word |= cp[level];
          if constexpr (!is_tree) {
            if (cp.code_length > level + 1) {
              if (cp[level]) { text1.emplace_back(local_text[cur_pos + i]); }
              else { text0.emplace_back(local_text[cur_pos + i]); }
            }
          }
        }
        if constexpr (!is_tree) {
          _zeros[level] += __builtin_popcountll(~word);
        }
        bv[level][cur_pos >> 6] = word;
      }
      //std::cout << "\n";

      // if there are remaining odd bits...
      if (local_text.size() & 63ULL) {
        uint64_t word = 0ULL;
        for (uint32_t i = 0; i < local_text.size() - cur_pos; ++i) {
          const code_pair cp = codes.encode_symbol(local_text[cur_pos + i]);
          word <<= 1;
          word |= cp[level];
          if constexpr (!is_tree) {
            if (cp.code_length > level + 1) {
              if (cp[level]) { text1.emplace_back(local_text[cur_pos + i]); }
              else { text0.emplace_back(local_text[cur_pos + i]); }
            }
          }
        }
        word <<= (64 - (local_text.size() & 63ULL));
        bv[level][local_text.size() >> 6] = word;
        if constexpr (!is_tree) {
          _zeros[level] +=
            (local_text.size() & 63ULL) - __builtin_popcountll(word);
        }
      }

      if constexpr (is_tree) {
        std::vector<std::vector<AlphabetType>> buckets(1ULL << (level + 1));

        for (uint64_t i = 0; i < local_text.size(); ++i) {
          const AlphabetType cur_symbol = local_text[i];
          if (codes.code_length(cur_symbol) > level + 1) {
            // prefix gets the prefix bits of the codeword
            buckets[codes.encode_symbol(cur_symbol).prefix(level + 1)]
              .emplace_back(cur_symbol);
          }
        }

        /*
        print_list(std::cout, buckets, true, [&](auto& out, auto& val) {
          print_list(out, val);
        }) << "\n";
        */

        for (uint64_t i = 1; i < buckets.size(); ++i) {
          std::move(buckets[i].begin(), buckets[i].end(),
            std::back_inserter(buckets[0]));
        }
        local_text.swap(buckets[0]);
      } else /*if constexpr (!is_tree)*/ {
        std::move(text1.begin(), text1.end(), std::back_inserter(text0));
        local_text = std::move(text0);
      }
      //print_list(std::cout, local_text) << "\n";
    }
    if constexpr (is_tree) {
      return wavelet_structure_tree_huffman<AlphabetType>(
        std::move(_bv), std::move(codes));
    } else /*if constexpr (!is_tree)*/ {
      return wavelet_structure_matrix_huffman<AlphabetType>(
        std::move(_bv), std::move(_zeros), std::move(codes));
    }
  }
}; // class wx_huff_naive

/******************************************************************************/
