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

#include "pc.hpp"

template <typename AlphabetType, typename SizeType>
class wt_dd_pc {

public:
    wt_dd_pc(const std::vector<AlphabetType>& global_text,
             const SizeType size,
             const SizeType levels) {
        if(global_text.size() == 0) { return; }

        const SizeType shards = omp_get_max_threads();

        // Do all bulk allocations in the same thread:
        // TODO: flatten vectors where possible, to reduce indirection

        const auto rho = rho_identity(levels);
        auto ctxs = std::vector<KeepLevel<SizeType, false, decltype(rho)>>(shards);

        for (size_t shard = 0; shard < shards; shard++) {
            const SizeType local_size = (size / shards) +
                ((shard < size % shards) ? 1 : 0);

            // TODO: store size in ctx so that later can reload it
            ctxs[shard] = KeepLevel<SizeType, false, decltype(rho)>(local_size, levels, rho);
        }

        #pragma omp parallel
        {
            const SizeType omp_rank = omp_get_thread_num();
            const SizeType omp_size = omp_get_num_threads();
            assert(omp_size == shards);

            const SizeType local_size = (size / omp_size) +
                ((omp_rank < size % omp_size) ? 1 : 0);
            const SizeType offset = (omp_rank * (size / omp_size)) +
                std::min<SizeType>(omp_rank, size % omp_size);

            const AlphabetType* text = global_text.data() + offset;

            pc(text, local_size, levels, ctxs[omp_rank]);
        }

        // TODO: Move and drop unneeded ctx stuff better than this

        auto glob_bv = std::vector<Bvs<SizeType>>(
            shards);

        auto glob_hist = std::vector<std::vector<std::vector<SizeType>>>(
            shards, std::vector<std::vector<SizeType>>(
                levels + 1, std::vector<SizeType>((1 << levels), 0)));

        for(size_t shard = 0; shard < shards; shard++) {
            glob_bv[shard] = std::move(ctxs[shard].bv());
            for(size_t level = 0; level < (levels + 1); level++) {
                auto hist_size = ctxs[shard].hist_size(level);
                glob_hist[shard][level] = std::vector<SizeType>(hist_size, 0);
                for(size_t i = 0; i < hist_size; i++) {
                    glob_hist[shard][level][i] = ctxs[shard].hist(level, i);
                }
            }
        }

        _bv = merge_bvs<SizeType>(size, levels, shards, glob_hist, glob_bv, rho);
    }

    auto get_bv() const {
        return _bv.vec();
    }

private:
    Bvs<SizeType> _bv;
}; // class wt_dd_pc

#endif // WT_DOMAIN_DECOMPOSITION_PREFIX_COUNTING

/******************************************************************************/
