/*******************************************************************************
 * include/wt_dd_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WT_DOMAIN_DECOMPOSITION_PREFIX_COUNTING
#define WT_DOMAIN_DECOMPOSITION_PREFIX_COUNTING

#include <algorithm>
#include <chrono>
#include <cstring>
#include <omp.h>
#include <vector>

#include "common.hpp"
#include "debug.hpp"
#include "merge.hpp"

template <typename AlphabetType, typename SizeType>
class wt_dd_pc {

public:
    wt_dd_pc(const std::vector<AlphabetType>& text,
             const SizeType size,
             const SizeType levels) {

        if(text.size() == 0) { return; }

        const SizeType shards = omp_get_max_threads();

        // Do all bulk allocations in the same thread:
        // TODO: flatten vectors where possible, to reduce indirection

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

        #pragma omp parallel
        {
            const size_t omp_rank = omp_get_thread_num();
            const size_t omp_size = omp_get_num_threads();

            assert(omp_size == shards);

            const SizeType local_size = (size / omp_size) +
                ((omp_rank < size % omp_size) ? 1 : 0);
            const SizeType offset = (omp_rank * (size / omp_size)) +
                std::min<SizeType>(omp_rank, size % omp_size);

            SizeType cur_max_char = (1 << levels);

            auto& borders = glob_borders[omp_rank];
            auto& hist = glob_hist[omp_rank];
            auto& bv = glob_bv[omp_rank];

            // While initializing the histogram, we also compute the fist level
            SizeType cur_pos = 0;
            for (; cur_pos + 64 <= local_size; cur_pos += 64) {
                uint64_t word = 0ULL;
                for (SizeType i = offset; i < 64 + offset; ++i) {
                    ++hist[levels][text[cur_pos + i]];
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                bv.vec()[0][cur_pos >> 6] = word;
            }
            if (local_size & 63ULL) {
                uint64_t word = 0ULL;
                for (SizeType i = offset; i < local_size - cur_pos + offset; ++i) {
                    ++hist[levels][text[cur_pos + i]];
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                word <<= (64 - (local_size & 63ULL));
                bv.vec()[0][local_size >> 6] = word;
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
                    borders[i] = borders[i - 1] + hist[level][i - 1];
                }

                // Now we insert the bits with respect to their bit prefixes
                for (SizeType i = offset; i < local_size + offset; ++i) {
                    const SizeType pos = borders[text[i] >> prefix_shift]++;
                    bv.vec()[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
                        << (63ULL - (pos & 63ULL)));
                }
            }

            hist[0].resize(1);
            if (hist.size() > 1) {
                hist[0][0] = hist[1][0] + hist[1][1];
            }
        }

        drop_me(std::move(glob_borders));

        // TODO: Add abstraction for allocating the bitvector (no more bare vector of pointers)

        const auto rho = rho_identity(levels);

        auto bv = merge_bvs<SizeType>(size, levels, shards, glob_hist, glob_bv, rho);
        _bv = std::move(bv.vec());

    }

    auto get_bv() const {
        return _bv;
    }

private:
    std::vector<uint64_t*> _bv;
}; // class wt_dd_pc

#endif // WT_DOMAIN_DECOMPOSITION_PREFIX_COUNTING

/******************************************************************************/