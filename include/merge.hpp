#pragma once

#include "common.hpp"
#include "debug.hpp"
#include <cassert>

template<typename SizeType, typename Rho>
inline auto merge_bvs(SizeType size,
                      SizeType levels,
                      SizeType shards,
                      const std::vector<std::vector<std::vector<SizeType>>>& glob_hist,
                      const std::vector<Bvs<SizeType>>& glob_bv,
                      const Rho& rho) -> Bvs<SizeType>
{
    assert(shards == glob_bv.size());
    assert(shards == glob_hist.size());

    auto offsets = std::vector<SizeType>(shards);
    auto local_offsets = std::vector<std::vector<std::vector<SizeType>>>(
        levels, std::vector<std::vector<SizeType>>(
            shards, std::vector<SizeType>(shards + 1)));
    auto word_offsets = std::vector<std::vector<SizeType>>(
        levels, std::vector<SizeType>(shards + 1));
    auto block_seq_offsets = std::vector<std::vector<SizeType>>(
        levels, std::vector<SizeType>(shards + 1));
    auto glob_cursors = std::vector<std::vector<std::vector<size_t>>>(
        shards, std::vector<std::vector<size_t>>(
            levels, std::vector<size_t>(shards)
        ));

    for (size_t rank = 1; rank < shards; rank++) {
        const size_t omp_rank = rank;
        const size_t omp_size = shards;

        const SizeType offset = (omp_rank * (word_size(size) / omp_size)) +
            std::min<SizeType>(omp_rank, word_size(size) % omp_size);

        offsets[rank - 1] = offset * 64ull;
    }
    offsets[shards - 1] = word_size(size) * 64ull;

    // TODO: fix
    for(size_t level = 0; level < levels; level++) {
        const size_t br_size = 1ull << level;
        size_t j = 0;
        size_t oi = 0;

        for(size_t i = 0; i < br_size * shards; i++) {
            const auto bit_i = i / shards;
            const auto shard = i % shards;

            const auto h = glob_hist[shard][level];
            const auto block_size = h[rho(level, bit_i)];

            j += block_size;
            local_offsets[level][shard][oi + 1] += block_size;

            if (j <= offsets[oi]) {
                // ...
            } else {
                size_t right = j - offsets[oi];
                size_t left = block_size - right;

                j -= block_size;
                local_offsets[level][shard][oi + 1] -= block_size;

                j += left;
                local_offsets[level][shard][oi + 1] += left;

                // TODO: rename word_offset to something like
                // "offset in block"
                word_offsets[level][oi + 1] = left;
                block_seq_offsets[level][oi + 1] = i;

                if (oi + 2 < shards + 1) {
                    for(size_t s = 0; s < shards; s++) {
                        local_offsets[level][s][oi + 2] = local_offsets[level][s][oi + 1];
                    }
                }
                oi++;

                j += right;
                local_offsets[level][shard][oi + 1] += right;
            }
        }
    }

    auto r = Bvs<SizeType>(size, levels);
    auto& _bv = r.vec();

    #pragma omp parallel
    {
        const size_t merge_shard = omp_get_thread_num();

        assert(size_t(omp_get_num_threads()) == shards);

        const auto target_right = std::min(offsets[merge_shard], size);
        const auto target_left = std::min((merge_shard > 0 ? offsets[merge_shard - 1] : 0), target_right);

        auto& cursors = glob_cursors[merge_shard];

        for (size_t level = 0; level < levels; level++) {
            for(size_t read_shard = 0; read_shard < shards; read_shard++) {
                cursors[level][read_shard] = local_offsets[level][read_shard][merge_shard];
            }
        }

        for (size_t level = 0; level < levels; level++) {
            auto seq_i = block_seq_offsets[level][merge_shard];

            size_t j = target_left;
            size_t init_offset = word_offsets[level][merge_shard];

            while (true) {
                const auto i = seq_i / shards;
                const auto shard = seq_i % shards;

                const auto& h = glob_hist[shard][level];
                const auto& local_bv = glob_bv[shard].vec()[level];

                auto block_size = h[rho(level, i)] - init_offset;
                init_offset = 0; // TODO: remove this by doing a initial pass
                auto& local_cursor = cursors[level][shard];

                // TODO: copy over whole block
                while(block_size != 0) {
                    if (j >= target_right) {
                        break;
                    }

                    block_size--;
                    const auto src_pos = local_cursor++;
                    const auto pos = j++;
                    const bool bit = bit_at(local_bv, src_pos);

                    _bv[level][pos >> 6] |= (uint64_t(bit) << (63ULL - (pos & 63ULL)));
                }

                if (j >= target_right) {
                    break;
                }

                seq_i++;
            }
        }
    }

    return r;
}
