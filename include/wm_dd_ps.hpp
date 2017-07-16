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
#include "merge.hpp"

template <typename AlphabetType, typename SizeType>
class wm_dd_ps {

public:
    wm_dd_ps(const std::vector<AlphabetType>& global_text,
             const SizeType size,
             const SizeType levels):
        _zeros(levels, 0)
    {

        if(global_text.size() == 0) { return; }

        const SizeType shards = omp_get_max_threads();

        // Do all bulk allocations in the same thread:
        // TODO: flatten vectors where possible, to reduce indirection

        auto glob_zeros = std::vector<std::vector<SizeType>>(
            shards, std::vector<SizeType>(levels));
        auto glob_bv = std::vector<Bvs<SizeType>>(
            shards);
        auto glob_hist = std::vector<std::vector<std::vector<SizeType>>>(
            shards, std::vector<std::vector<SizeType>>(
                levels + 1, std::vector<SizeType>((1 << levels), 0)));
        auto glob_borders = std::vector<std::vector<SizeType>>(
            shards, std::vector<SizeType>((1 << levels), 0));

        for (size_t shard = 0; shard < shards; shard++) {
            const SizeType local_size = (size / shards) +
                ((shard < size % shards) ? 1 : 0);

            glob_bv[shard] = Bvs<SizeType>(local_size, levels);
        }

        std::vector<AlphabetType> global_sorted_text(size);

        #pragma omp parallel
        {
            const SizeType omp_rank = omp_get_thread_num();
            const SizeType omp_size = omp_get_num_threads();

            assert(omp_size == shards);

            const SizeType local_size = (size / omp_size) +
                ( (omp_rank < size % omp_size) ? 1 : 0);
            const SizeType offset = (omp_rank * (size / omp_size)) +
                std::min<SizeType>(omp_rank, size % omp_size);

            const AlphabetType* text = global_text.data() + offset;
            AlphabetType* sorted_text = global_sorted_text.data() + offset;

            SizeType cur_max_char = (1 << levels);
            std::vector<SizeType> bit_reverse = BitReverse<SizeType>(levels - 1);

            auto& zeros = glob_zeros[omp_rank];
            auto& borders = glob_borders[omp_rank];
            auto& hist = glob_hist[omp_rank];
            auto& bv = glob_bv[omp_rank];

            // While initializing the histogram, we also compute the fist level
            SizeType cur_pos = 0;
            for (; cur_pos + 64 <= local_size; cur_pos += 64) {
                uint64_t word = 0ULL;
                for (SizeType i = 0; i < 64; ++i) {
                    ++hist[levels][text[cur_pos + i]];
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                bv.vec()[0][cur_pos >> 6] = word;
            }
            if (local_size & 63ULL) {
                uint64_t word = 0ULL;
                for (SizeType i = 0; i < local_size - cur_pos; ++i) {
                    ++hist[levels][text[cur_pos + i]];
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                word <<= (64 - (local_size & 63ULL));
                bv.vec()[0][local_size >> 6] = word;
            }

            // The number of 0s at the last level is the number of "even" characters
            for (SizeType i = 0; i < cur_max_char; i += 2) {
                zeros[levels - 1] += hist[levels][i];
            }

            for (SizeType level = levels - 1; level > 0; --level) {
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
                    bv.vec()[level][cur_pos >> 6] = word;
                }
                if (local_size & 63ULL) {
                    uint64_t word = 0ULL;
                    for (SizeType i = 0; i < local_size - cur_pos; ++i) {
                        word <<= 1;
                        word |= ((sorted_text[cur_pos + i] >> ((levels - 1) - level)) & 1ULL);
                    }
                    word <<= (64 - (local_size & 63ULL));
                    bv.vec()[level][local_size >> 6] = word;
                }
            }

            hist[0].resize(1);
            if (hist.size() > 1) {
                hist[0][0] = hist[1][0] + hist[1][1];
            }
        }

        // TODO: Add abstraction for allocating the bitvector (no more bare vector of pointers)

        const auto rho = rho_bit_reverse(levels);

        auto bv = merge_bvs<SizeType>(size, levels, shards, glob_hist, glob_bv, rho);
        _bv = std::move(bv.vec());

        #pragma omp parallel for
        for(size_t level = 0; level < levels; level++) {
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
}; // class wm_dd_ps

#endif // WM_DOMAIN_DECOMPOSITION_PREFIX_SORTING

/******************************************************************************/
