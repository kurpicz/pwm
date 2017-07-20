#pragma once

#include "common.hpp"
#include "debug.hpp"
#include <cassert>
#include <climits>
#include <omp.h>

template<typename SizeType, typename WordType>
void copy_bits(WordType* const dst,
               WordType const* const src,
               SizeType& dst_off_ref,
               SizeType& src_off_ref,
               SizeType const block_size)
{
    if (block_size == 0) return;

    WordType constexpr BITS = (sizeof(WordType) * CHAR_BIT);
    WordType constexpr MOD_MASK = BITS - 1;
    WordType constexpr SHIFT = log2(MOD_MASK);

    SizeType dst_off = dst_off_ref;
    SizeType src_off = src_off_ref;
    auto const dst_off_end = dst_off + block_size;

    // NB: Check if source and target block are aligned the same way
    // This is needed because the non-aligned code path would
    // trigger undefined behavior in the aligned case.
    if ((dst_off & MOD_MASK) == (src_off & MOD_MASK)) {
        // Copy unaligned leading bits
        {
            auto& word = dst[dst_off >> SHIFT];

            while ((dst_off & MOD_MASK) != 0 && dst_off != dst_off_end) {
                bool const bit = bit_at<WordType>(src, src_off++);

                word |= (WordType(bit) << (MOD_MASK - (dst_off++ & MOD_MASK)));
            }
        }

        // Copy the the bulk in-between word-wise
        {
            auto const words = (dst_off_end - dst_off) >> SHIFT;

            WordType*       ds = dst + (dst_off >> SHIFT);
            WordType const* sr = src + (src_off >> SHIFT);

            WordType const* const ds_end = ds + words;

            while (ds != ds_end) {
                *ds++ = *sr++;
            }

            dst_off += words * BITS;
            src_off += words * BITS;
        }

        // Copy unaligned trailing bits
        {
            auto& word = dst[dst_off >> SHIFT];

            while (dst_off != dst_off_end) {
                bool const bit = bit_at<WordType>(src, src_off++);

                word |= (WordType(bit) << (MOD_MASK - (dst_off++ & MOD_MASK)));
            }
        }
    } else {
        // Copy unaligned leading bits
        {
            auto& word = dst[dst_off >> SHIFT];

            while ((dst_off & MOD_MASK) != 0 && dst_off != dst_off_end) {
                bool const bit = bit_at<WordType>(src, src_off++);

                word |= (WordType(bit) << (MOD_MASK - (dst_off++ & MOD_MASK)));
            }
        }

        // Copy the the bulk in-between word-wise
        {
            auto const words = (dst_off_end - dst_off) >> SHIFT;

            WordType const src_shift_a = src_off & MOD_MASK;
            WordType const src_shift_b = BITS - src_shift_a;

            WordType*       ds = dst + (dst_off >> SHIFT);
            WordType const* sr = src + (src_off >> SHIFT);
            WordType const* const ds_end = ds + words;

            while (ds != ds_end) {
                *ds++ = (*sr << src_shift_a) | (*(sr+1) >> src_shift_b);
                sr++;
            }

            dst_off += words * BITS;
            src_off += words * BITS;
        }

        // Copy unaligned trailing bits
        {
            auto& word = dst[dst_off >> SHIFT];

            while (dst_off != dst_off_end) {
                bool const bit = bit_at<WordType>(src, src_off++);

                word |= (WordType(bit) << (MOD_MASK - (dst_off++ & MOD_MASK)));
            }
        }
    }

    dst_off_ref += block_size;
    src_off_ref += block_size;
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

    // Initialize bulk data structures centrally

    struct MergeCtx {
        SizeType offset;
        std::vector<std::vector<SizeType>> read_offsets;
        std::vector<SizeType> offsets_in_first_word;
        std::vector<SizeType> first_read_block;
    };

    auto ctxs = std::vector<MergeCtx>(shards, MergeCtx {
        0,
        std::vector<std::vector<SizeType>>(
            levels, std::vector<SizeType>(shards)),
        std::vector<SizeType>(
            levels),
        std::vector<SizeType>(
            levels),
    });

    // Calculate bit offset per merge shard (thread)
    for (size_t rank = 1; rank < shards; rank++) {
        const size_t omp_rank = rank;
        const size_t omp_size = shards;

        const SizeType offset = (omp_rank * (word_size(size) / omp_size)) +
            std::min<SizeType>(omp_rank, word_size(size) % omp_size);

        ctxs[rank - 1].offset = offset * 64ull;
    }
    ctxs[shards - 1].offset = word_size(size) * 64ull;

    for(size_t level = 0; level < levels; level++) {
        const size_t br_size = 1ull << level;

        size_t write_offset = 0; // bit offset in destination bv
        size_t merge_shard = 0; // index of merge thread

        for(size_t i = 0; i < br_size * shards; i++) {
            const auto bit_i = i / shards;
            const auto read_shard = i % shards;

            auto read_offset = [&level, &ctxs](auto merge_shard, auto read_shard) -> SizeType& {
                return ctxs[merge_shard].read_offsets[level][read_shard];
            };

            const auto h = glob_hist[read_shard][level];

            auto block_size = h[rho(level, bit_i)];

            write_offset += block_size;
            read_offset(merge_shard + 1, read_shard) += block_size;

            // If we passed the current right border, split up the block
            if (write_offset > ctxs[merge_shard].offset) {
                // Take back the last step
                write_offset -= block_size;
                read_offset(merge_shard + 1, read_shard) -= block_size;

                SizeType offset_in_first_word = 0;
                do {
                    // Split up the block like this:
                    //       [    left_block_size    |        right_block_size        ]
                    //       ^                       ^                                ^
                    // (write_offset)    (ctxs[merge_shard].offset)      (write_offset + block_size)

                    auto const left_block_size = ctxs[merge_shard].offset - write_offset;

                    write_offset += left_block_size;
                    read_offset(merge_shard + 1, read_shard) += left_block_size;

                    offset_in_first_word += left_block_size;
                    ctxs[merge_shard + 1].offsets_in_first_word[level] = offset_in_first_word;
                    ctxs[merge_shard + 1].first_read_block[level] = i;

                    if (merge_shard + 2 < shards) {
                        for(size_t s = 0; s < shards; s++) {
                            read_offset(merge_shard + 2, s) = read_offset(merge_shard + 1, s);
                        }
                    }

                    merge_shard++;

                    // Once we have calculated the offsets for all merge threads,
                    // break out of the whole nested loop
                    if (merge_shard + 1 == shards) {
                        goto triple_loop_exit;
                    }

                    // Iterate on remaining block, because one block might
                    // span multiple threads
                    block_size -= left_block_size;
                } while ((write_offset + block_size) > ctxs[merge_shard].offset);

                // Process remainder of block
                write_offset += block_size;
                read_offset(merge_shard + 1, read_shard) += block_size;

                assert(write_offset <= ctxs[merge_shard].offset);
            }
        }
        triple_loop_exit:; // we are done
    }

    auto r = Bvs<SizeType>(size, levels);
    auto& _bv = r.vec();

    #pragma omp parallel
    //for(size_t omp_rank = 0; omp_rank < shards; omp_rank++)
    {
        assert(size_t(omp_get_num_threads()) == shards);
        const size_t omp_rank = omp_get_thread_num();
        const size_t merge_shard = omp_rank;

        auto& ctx = ctxs[merge_shard];

        const auto target_right = std::min(ctx.offset, size);
        const auto target_left = std::min((merge_shard > 0 ? ctxs[merge_shard - 1].offset : 0), target_right);

        for (size_t level = 0; level < levels; level++) {
            auto seq_i = ctx.first_read_block[level];

            SizeType write_offset = target_left;
            size_t init_offset = ctx.offsets_in_first_word[level];

            while (write_offset < target_right) {
                const auto i = seq_i / shards;
                const auto shard = seq_i % shards;
                seq_i++;

                const auto& h = glob_hist[shard][level];
                const auto& local_bv = glob_bv[shard].vec()[level];

                auto const block_size = std::min<SizeType>(
                    target_right - write_offset,
                    h[rho(level, i)] - init_offset
                );
                init_offset = 0; // TODO: remove this by doing a initial pass

                auto& local_cursor = ctx.read_offsets[level][shard];

                copy_bits<SizeType, uint64_t>(
                    _bv[level],
                    local_bv,
                    write_offset,
                    local_cursor,
                    block_size
                );
            }
        }
    }

    return r;
}
