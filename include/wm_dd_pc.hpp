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
#include <cstring>
#include <omp.h>
#include <vector>

#include "common.hpp"
#include "merge.hpp"
#include "pc.hpp"

template <typename AlphabetType, typename SizeType>
class wm_dd_pc {

public:
    wm_dd_pc(const std::vector<AlphabetType>& global_text,
             const SizeType size,
             const SizeType levels):
        _zeros(levels, 0)
    {

        if(global_text.size() == 0) { return; }

        const SizeType shards = omp_get_max_threads();

        // Do all bulk allocations in the same thread:
        // TODO: flatten vectors where possible, to reduce indirection

        const auto rho = rho_bit_reverse(levels);
        auto ctxs = std::vector<KeepLevel<SizeType, true, decltype(rho)>>(shards);

        for (size_t shard = 0; shard < shards; shard++) {
            const SizeType local_size = (size / shards) +
                ((shard < size % shards) ? 1 : 0);

            // TODO: store size in ctx so that later can reload it
            ctxs[shard] = KeepLevel<SizeType, true, decltype(rho)>(local_size, levels, rho);
        }

        #pragma omp parallel
        {
            const SizeType omp_rank = omp_get_thread_num();
            const SizeType omp_size = omp_get_num_threads();
            assert(omp_size == shards);

            auto& ctx = ctxs[omp_rank];

            const SizeType local_size = (size / omp_size) +
                ((omp_rank < size % omp_size) ? 1 : 0);
            const SizeType offset = (omp_rank * (size / omp_size)) +
                std::min<SizeType>(omp_rank, size % omp_size);

            const AlphabetType* text = global_text.data() + offset;

            SizeType cur_max_char = (1 << levels);

            auto& zeros = ctx.zeros();
            auto& borders = ctx.borders();
            auto& bv = ctx.bv().vec();

            // While initializing the histogram, we also compute the fist level
            SizeType cur_pos = 0;
            for (; cur_pos + 64 <= local_size; cur_pos += 64) {
                uint64_t word = 0ULL;
                for (SizeType i = 0; i < 64; ++i) {
                    ++ctx.hist(levels, text[cur_pos + i]);
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                bv[0][cur_pos >> 6] = word;
            }
            if (local_size & 63ULL) {
                uint64_t word = 0ULL;
                for (SizeType i = 0; i < local_size - cur_pos; ++i) {
                    ++ctx.hist(levels, text[cur_pos + i]);
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                word <<= (64 - (local_size & 63ULL));
                bv[0][local_size >> 6] = word;
            }

            // The number of 0s at the last level is the number of "even" characters
            if (ctx.calc_zeros()) {
                for (SizeType i = 0; i < cur_max_char; i += 2) {
                    zeros[levels - 1] += ctx.hist(levels, i);
                }
            }

            // Now we compute the WM bottom-up, i.e., the last level first
            for (SizeType level = levels - 1; level > 0; --level) {
                const SizeType prefix_shift = (levels - level);
                const SizeType cur_bit_shift = prefix_shift - 1;

                // Update the maximum value of a feasible a bit prefix and update the
                // histogram of the bit prefixes
                cur_max_char >>= 1;

                for (SizeType i = 0; i < cur_max_char; ++i) {
                    ctx.hist(level, i)
                        = ctx.hist(level + 1, i << 1)
                        + ctx.hist(level + 1, (i << 1) + 1);
                }

                // Compute the starting positions of characters with respect to their
                // bit prefixes and the bit-reversal permutation
                borders[0] = 0;
                for (SizeType i = 1; i < cur_max_char; ++i) {
                    auto const prev_rho = ctx.rho(level, i - 1);

                    borders[ctx.rho(level, i)] = borders[prev_rho] + ctx.hist(level, prev_rho);

                    if (ctx.calc_rho())  {
                        ctx.set_rho(level - 1, i - 1, prev_rho >> 1);
                    }
                }

                // The number of 0s is the position of the first 1 in the previous level
                if (ctx.calc_zeros()) {
                    zeros[level - 1] = borders[1];
                }

                // Now we insert the bits with respect to their bit prefixes
                for (SizeType i = 0; i < local_size; ++i) {
                    const SizeType pos = borders[text[i] >> prefix_shift]++;
                    bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
                        << (63ULL - (pos & 63ULL)));
                }
            }

            if (levels > 1) { // TODO check condition
                ctx.hist(0, 0) = ctx.hist(1, 0) + ctx.hist(1, 1);
            }
        }

        auto glob_bv = std::vector<Bvs<SizeType>>(
            shards);

        auto glob_zeros = std::vector<std::vector<SizeType>>(
            shards, std::vector<SizeType>(levels));

        auto glob_hist = std::vector<std::vector<std::vector<SizeType>>>(
            shards, std::vector<std::vector<SizeType>>(
                levels + 1, std::vector<SizeType>((1 << levels), 0)));

        for(size_t shard = 0; shard < shards; shard++) {
            glob_bv[shard] = std::move(ctxs[shard].bv());
            glob_zeros[shard] = std::move(ctxs[shard].zeros());
            for(size_t level = 0; level < (levels + 1); level++) {
                auto hist_size = ctxs[shard].hist_size(level);
                glob_hist[shard][level] = std::vector<SizeType>(hist_size, 0);
                for(size_t i = 0; i < hist_size; i++) {
                    glob_hist[shard][level][i] = ctxs[shard].hist(level, i);
                }
            }
        }

        // TODO: Add abstraction for allocating the bitvector (no more bare vector of pointers)

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
}; // class wm_dd_pc

#endif // WM_DOMAIN_DECOMPOSITION_PREFIX_COUNTING

/******************************************************************************/
