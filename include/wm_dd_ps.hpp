/*******************************************************************************
 * include/wm_dd_ps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_DOMAIN_DECOMPOSITION_PREFIX_SORTING
#define WM_DOMAIN_DECOMPOSITION_PREFIX_SORTING

#include <algorithm>
#include <cstring>
#include <omp.h>
#include <vector>

#include "common.hpp"

template <typename AlphabetType, typename SizeType>
class wm_dd_ps {

public:
    wm_dd_ps(const std::vector<AlphabetType>& global_text,
             const SizeType size,
             const SizeType levels):
        _bv(levels), _zeros(levels, 0)
    {

        if(global_text.size() == 0) { return; }

        std::vector<AlphabetType> global_sorted_text(size);

        #pragma omp parallel
        {
            const auto omp_rank = omp_get_thread_num();
            const auto omp_size = omp_get_num_threads();

            const SizeType local_size = (size / omp_size) +
                ( (omp_rank < size % omp_size) ? 1 : 0);
            const SizeType offset = (omp_rank * (size / omp_size)) +
                std::min(static_cast<uint32_t>(omp_rank), size % omp_size);

            const AlphabetType* text = global_text.data() + offset;
            AlphabetType* sorted_text = global_sorted_text.data() + offset;

            SizeType cur_max_char = (1 << levels);
            std::vector<SizeType> util_vec(cur_max_char << 1, 0);
            std::vector<SizeType> bit_reverse = BitReverse<SizeType>(levels - 1);
            SizeType* hist = util_vec.data();
            SizeType* borders = util_vec.data() + cur_max_char;

            std::vector<uint64_t*> tmp_bv(levels);
            std::vector<size_t> tmp_zeros(levels, 0);
            tmp_bv[0] = new uint64_t[((local_size + 63ULL) >> 6) * levels];
            // memset is ok (all to 0)
            memset(tmp_bv[0], 0, (((local_size + 63ULL) >> 6) * sizeof(uint64_t)) * levels);
            for (SizeType level = 1; level < levels; ++level) {
                tmp_bv[level] = tmp_bv[level - 1] + ((local_size + 63ULL) >> 6);
            }

            // While initializing the histogram, we also compute the fist level
            SizeType cur_pos = 0;
            for (; cur_pos + 64 <= local_size; cur_pos += 64) {
                uint64_t word = 0ULL;
                for (SizeType i = 0; i < 0 + 64; ++i) {
                    ++hist[text[cur_pos + i]];
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                (tmp_bv[0])[cur_pos >> 6] = word;
            }
            if (local_size & 63ULL) {
                uint64_t word = 0ULL;
                for (SizeType i = 0; i < local_size - cur_pos; ++i) {
                    ++hist[text[cur_pos + i]];
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                word <<= (64 - (local_size & 63ULL));
                (tmp_bv[0])[local_size >> 6] = word;
            }

            // The number of 0s at the last level is the number of "even" characters
            for (SizeType i = 0; i < cur_max_char; i += 2) {
                tmp_zeros[levels - 1] += hist[i];
            }

            for (SizeType level = levels - 1; level > 0; --level) {
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
                tmp_zeros[level - 1] = borders[1];

                // Now we sort the text utilizing counting sort and the starting positions
                // that we have computed before
                for (SizeType i = 0; i < local_size; ++i) {
                    const AlphabetType cur_char = text[i];
                    sorted_text[borders[cur_char >> (levels - level)]++] = cur_char;
                }

                // Since we have sorted the text, we can simply scan it from left to right
                // and for the character at position $i$ we set the $i$-th bit in the
                // bit vector accordingly
                for (cur_pos = 0; cur_pos + 63 < local_size; cur_pos += 64) {
                    uint64_t word = 0ULL;
                    for (SizeType i = 0; i < 64; ++i) {
                        word <<= 1;
                        word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
                    }
                    tmp_bv[level][cur_pos >> 6] = word;
                }
                if (local_size & 63ULL) {
                    uint64_t word = 0ULL;
                    for (SizeType i = 0; i < local_size - cur_pos; ++i) {
                        word <<= 1;
                        word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
                    }
                    word <<= (64 - (local_size & 63ULL));
                    tmp_bv[level][local_size >> 6] = word;
                }
            }
        }
    }

    auto get_bv_and_zeros() const {
        return std::make_pair(_bv, _zeros);
    }

private:
    std::vector<uint64_t*> _bv;
    std::vector<SizeType> _zeros;
}; // class wm_dd_ps

#endif // WM_DOMAIN_DECOMPOSITION_PREFIX_SORTING

/******************************************************************************/
