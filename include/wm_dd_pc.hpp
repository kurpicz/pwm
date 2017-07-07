/*******************************************************************************
 * include/wm_dd_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_DOMAIN_DECOMPOSITION_PREFIX_COUNTING
#define WM_DOMAIN_DECOMPOSITION_PREFIX_COUNTING

#include <algorithm>
#include <chrono>
#include <cstring>
#include <omp.h>
#include <vector>

#include "common.hpp"

template <typename AlphabetType, typename SizeType>
class wm_dd_pc {

public:
  wm_dd_pc(const std::vector<AlphabetType>& text, const SizeType size,
    const SizeType levels) : _bv(levels), _zeros(levels, 0) {

    if(text.size() == 0) { return; }

    std::chrono::system_clock::time_point first, second, third, fourth;

    #pragma omp parallel
    {
      const size_t omp_rank = omp_get_thread_num();
      const size_t omp_size = omp_get_num_threads();

      #pragma omp single
      first = std::chrono::system_clock::now();

      const SizeType local_size = (size / omp_size) +
        ((omp_rank < size % omp_size) ? 1 : 0);
      const SizeType offset = (omp_rank * (size / omp_size)) +
        std::min<SizeType>(omp_rank, size % omp_size);

      SizeType cur_max_char = (1 << levels);
      std::vector<SizeType> bit_reverse = BitReverse<SizeType>(levels - 1);
      std::vector<SizeType> hist(cur_max_char, 0);
      std::vector<SizeType> borders(cur_max_char, 0);

      std::vector<uint64_t*> tmp_bv(levels);

      tmp_bv[0] = new uint64_t[((local_size + 63ULL) >> 6) * levels];
      // memset is ok (all to 0)
      memset(tmp_bv[0], 0, (((local_size + 63ULL) >> 6) * sizeof(uint64_t)) * levels);
      for (SizeType level = 1; level < levels; ++level) {
        tmp_bv[level] = tmp_bv[level - 1] + ((local_size + 63ULL) >> 6);
      }

      #pragma omp single
      second = std::chrono::system_clock::now();
      // While initializing the histogram, we also compute the fist level
      SizeType cur_pos = 0;
      for (; cur_pos + 64 <= local_size; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (SizeType i = offset; i < 64 + offset; ++i) {
          ++hist[text[cur_pos + i]];
          word <<= 1;
          word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        tmp_bv[0][cur_pos >> 6] = word;
      }
      if (local_size & 63ULL) {
        uint64_t word = 0ULL;
        for (SizeType i = offset; i < local_size - cur_pos + offset; ++i) {
          ++hist[text[cur_pos + i]];
          word <<= 1;
          word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (local_size & 63ULL));
        tmp_bv[0][local_size >> 6] = word;
      }

      // The number of 0s at the last level is the number of "even" characters
      for (SizeType i = 0; i < cur_max_char; i += 2) {
        _zeros[levels - 1] += hist[i];
      }

      #pragma omp single
      third = std::chrono::system_clock::now();

      // Now we compute the WM bottom-up, i.e., the last level first
      for (SizeType level = levels - 1; level > 0; --level) {
        const SizeType prefix_shift = (levels - level);
        const SizeType cur_bit_shift = prefix_shift - 1;

        // Update the maximum value of a feasible a bit prefix and update the
        // histogram of the bit prefixes
        cur_max_char >>= 1;
        for (SizeType i = 0; i < cur_max_char; ++i) {
          hist[i] = hist[i << 1] + hist[(i << 1) + 1];
        }

        // Compute the starting positions of characters with respect to their
        // bit prefixes and the bit-reversal permutation
        borders[0] = 0;
        for (SizeType i = 1; i < cur_max_char; ++i) {
          borders[bit_reverse[i]] = borders[bit_reverse[i - 1]] +
            hist[bit_reverse[i - 1]];
          bit_reverse[i - 1] >>= 1;
        }
        // The number of 0s is the position of the first 1 in the previous level
        _zeros[level - 1] = borders[1];

        // Now we insert the bits with respect to their bit prefixes
        for (SizeType i = offset; i < local_size + offset; ++i) {
          const SizeType pos = borders[text[i] >> prefix_shift]++;
          tmp_bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
            << (63ULL - (pos & 63ULL)));
        }
      }
      #pragma omp single
      fourth = std::chrono::system_clock::now();
      #pragma omp single
      std::cout << "First: " << std::chrono::duration_cast<std::chrono::microseconds>(second - first).count() << std::endl
                << "Second: " <<  std::chrono::duration_cast<std::chrono::microseconds>(third - second).count() << std::endl
                << "Third: " << std::chrono::duration_cast<std::chrono::microseconds>(fourth - third).count() << std::endl;
    }
  }

  auto get_bv_and_zeros() const {
    return std::make_pair(_bv, _zeros);
  }

private:
  std::vector<uint64_t*> _bv;
  std::vector<SizeType> _zeros;
}; // class wm_dd_pc

#endif // WM_DOMAIN_DECOMPOSITION_PREFIX_COUNTING

/******************************************************************************/
