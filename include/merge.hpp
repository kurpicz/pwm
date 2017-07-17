#pragma once

#include "common.hpp"
#include "debug.hpp"
#include <cassert>

template<typename SizeType>
void copy_bits(uint64_t* dst,
               uint64_t const* src,
               SizeType& j,
               SizeType& local_cursor,
               SizeType const block_size) {

    auto const end_j = j + block_size;


    auto copy_bits = [&](auto end) {
        while(j != end) {
            const bool bit = bit_at(src, local_cursor++);
            const auto pos = j++;

            dst[pos >> 6] |= (uint64_t(bit) << (63ULL - (pos & 63ULL)));
        }
    };

    auto word_align = (64ull - (j % 64ull));

    if (block_size < word_align) {
        copy_bits(end_j);
    } else {

        // TODO: Use one word shift somehow
        copy_bits(j + word_align);

        // TODO: Reduce this to just "words"
        auto const word_end = (end_j >> 6 << 6);

        std::stringstream ss;

        ss << "j: "
            << (j - word_align)
            << " - "
            << j
            << " - "
            << word_end
            << " - "
            << end_j
            << "\n";
        std::cout << ss.str();

        auto const words = (word_end - j) >> 6;

        uint64_t* dst = &dst[j >> 6];
        uint64_t* const dst_end = dst + words;

        uint64_t const* src = &src[local_cursor >> 6];

        auto const shift = local_cursor % 64;

        while (dst != dst_end) {
            *dst++;
            *src << shift | *src++ >> (63ull - shift);
        }

        //j = word_end;
        //local_cursor += words * 64;

        copy_bits(word_end);
        copy_bits(end_j);
    }
}

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

            while (j < target_right) {
                const auto i = seq_i / shards;
                const auto shard = seq_i % shards;
                seq_i++;

                const auto& h = glob_hist[shard][level];
                const auto& local_bv = glob_bv[shard].vec()[level];

                auto const block_size = std::min(
                    target_right - j,
                    h[rho(level, i)] - init_offset
                );
                init_offset = 0; // TODO: remove this by doing a initial pass

                auto& local_cursor = cursors[level][shard];

                copy_bits(
                    _bv[level],
                    local_bv,
                    j,
                    local_cursor,
                    block_size
                );
            }
        }
    }

    return r;
}
