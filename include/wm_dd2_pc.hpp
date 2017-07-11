/*******************************************************************************
 * include/wm_dd2_pc.hpp
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
#include "debug.hpp"

template <typename AlphabetType, typename SizeType>
class wm_dd2_pc {

public:
    wm_dd2_pc(const std::vector<AlphabetType>& text,
             const SizeType size,
             const SizeType levels) : _bv(levels), _zeros(levels, 0) {

        if(text.size() == 0) { return; }

        // TODO: flatten vectors

        // This will be used for merging
        auto glob_zeros = std::vector<std::vector<SizeType>>(
            omp_get_max_threads(), std::vector<SizeType>(levels));
        auto glob_bv = std::vector<std::vector<uint64_t*>>(
            omp_get_max_threads()
        );
        auto glob_hist = std::vector<std::vector<std::vector<SizeType>>>(
            omp_get_max_threads(), std::vector<std::vector<SizeType>>(
                levels + 1, std::vector<SizeType>((1 << levels), 0)));

        // This just exists for central allocation
        auto glob_borders = std::vector<std::vector<SizeType>>(
            omp_get_max_threads(), std::vector<SizeType>((1 << levels), 0));

        #pragma omp parallel
        {
            const size_t omp_rank = omp_get_thread_num();
            const size_t omp_size = omp_get_num_threads();

            const SizeType local_size = (size / omp_size) +
                ((omp_rank < size % omp_size) ? 1 : 0);
            const SizeType offset = (omp_rank * (size / omp_size)) +
                std::min<SizeType>(omp_rank, size % omp_size);

            SizeType cur_max_char = (1 << levels);
            std::vector<SizeType> bit_reverse = BitReverse<SizeType>(levels - 1);

            auto& zeros = glob_zeros[omp_rank];
            auto& borders = glob_borders[omp_rank];
            auto& hist = glob_hist[omp_rank];
            auto& tmp_bv = glob_bv[omp_rank];

            tmp_bv = std::vector<uint64_t*>(levels);

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
                for (SizeType i = offset; i < 64 + offset; ++i) {
                    ++hist[levels][text[cur_pos + i]];
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                tmp_bv[0][cur_pos >> 6] = word;
            }
            if (local_size & 63ULL) {
                uint64_t word = 0ULL;
                for (SizeType i = offset; i < local_size - cur_pos + offset; ++i) {
                    ++hist[levels][text[cur_pos + i]];
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                word <<= (64 - (local_size & 63ULL));
                tmp_bv[0][local_size >> 6] = word;
            }

            // The number of 0s at the last level is the number of "even" characters
            for (SizeType i = 0; i < cur_max_char; i += 2) {
                zeros[levels - 1] += hist[levels][i];
            }

            // Now we compute the WM bottom-up, i.e., the last level first
            for (SizeType level = levels - 1; level > 0; --level) {
                const SizeType prefix_shift = (levels - level);
                const SizeType cur_bit_shift = prefix_shift - 1;

                // Update the maximum value of a feasible a bit prefix and update the
                // histogram of the bit prefixes
                cur_max_char >>= 1;
                for (SizeType i = 0; i < cur_max_char; ++i) {
                    hist[level][i] = hist[level + 1][i << 1] + hist[level + 1][(i << 1) + 1];
                }
                hist[level].resize(cur_max_char);

                // Compute the starting positions of characters with respect to their
                // bit prefixes and the bit-reversal permutation
                borders[0] = 0;
                for (SizeType i = 1; i < cur_max_char; ++i) {
                    borders[bit_reverse[i]] = borders[bit_reverse[i - 1]] +
                        hist[level][bit_reverse[i - 1]];
                    bit_reverse[i - 1] >>= 1;
                }

                // The number of 0s is the position of the first 1 in the previous level
                zeros[level - 1] = borders[1];

                // Now we insert the bits with respect to their bit prefixes
                for (SizeType i = offset; i < local_size + offset; ++i) {
                    const SizeType pos = borders[text[i] >> prefix_shift]++;
                    tmp_bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
                        << (63ULL - (pos & 63ULL)));
                }
            }

            hist[0].resize(1);
            if (hist.size() > 1) {
                hist[0][0] = hist[1][0] + hist[1][1];
            }
        }

        // TODO: Add abstraction for allocating the bitvector (no more bare vector of pointers)

        auto bit_reverse = std::vector<std::vector<SizeType>>(levels);
        bit_reverse[levels - 1] = BitReverse<SizeType>(levels - 1);
        for(size_t level = levels - 1; level > 0; level--) {
            bit_reverse[level - 1] = std::vector<SizeType>(bit_reverse[level].size() / 2);
            for(size_t i = 0; i < bit_reverse[level - 1].size(); i++) {
                bit_reverse[level - 1][i] = bit_reverse[level][i] >> 1;
            }
        }

        auto cursors = std::vector<std::vector<size_t>>(levels, std::vector<size_t>(glob_bv.size()));

        #pragma omp parallel for
        for(size_t level = 0; level < levels; level++) {
            const auto& br = bit_reverse[level];
            size_t j = 0;

            _bv[level] = new uint64_t[(size + 63ULL) >> 6];
            // memset is ok (all to 0)
            memset(_bv[level], 0, ((size + 63ULL) >> 6) * sizeof(uint64_t));

            for(size_t i = 0; i < br.size(); i++) {
                for(size_t shard = 0; shard < glob_bv.size(); shard++) {
                    auto& h = glob_hist[shard][level];
                    auto& block_size = h[br[i]];
                    const auto& local_bv = glob_bv[shard][level];
                    auto& local_cursor = cursors[level][shard];

                    // TODO: copy over whole block
                    while(block_size != 0) {
                        block_size--;
                        const auto src_pos = local_cursor++;
                        const auto pos = j++;
                        const bool bit = bit_at(local_bv, src_pos);

                        _bv[level][pos >> 6] |= (uint64_t(bit) << (63ULL - (pos & 63ULL)));
                    }
                }
            }
            for(size_t shard = 0; shard < glob_bv.size(); shard++) {
                _zeros[level] += glob_zeros[shard][level];
            }
        }
    }

    auto get_bv_and_zeros() const {
        return std::make_pair(_bv, _zeros);
    }

private:
    std::vector<uint64_t*> _bv;
    std::vector<SizeType> _zeros;
}; // class wm_dd2_pc

#endif // WM_DOMAIN_DECOMPOSITION_PREFIX_COUNTING

/******************************************************************************/
