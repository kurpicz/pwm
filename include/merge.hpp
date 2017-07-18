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
               SizeType const block_size,
               SizeType const words_size,
               SizeType const full_words_size
              )
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

            auto chk = [&]() {
                std::cout
                    << "dst ptr: " << (ds - dst)
                    << ", src ptr: " << (sr - src)
                    << ", len: " << full_words_size
                    << ", shifts: " << size_t(src_shift_a) << ", " << size_t(src_shift_b)
                    << "\n";
                assert((ds - dst) < full_words_size);
                assert((sr - src) < full_words_size);
                assert(((sr - src) + 1) < full_words_size);
            };

            while (ds != ds_end) {
                chk();

                WordType const vala = *sr;
                std::cout << "sp pre:  " << size_t(sr) << "\n";
                sr++;
                std::cout << "sp post: " << size_t(sr) << "\n";
                WordType const valb = *sr;

                *ds++ = (vala << src_shift_a) | (valb >> src_shift_b);

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

    //#pragma omp parallel
    for(size_t omp_rank = 0; omp_rank < shards; omp_rank++)
    {
        //assert(size_t(omp_get_num_threads()) == shards);
        //const size_t omp_rank = omp_get_thread_num();
        const size_t merge_shard = omp_rank;

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
                    block_size,
                    word_size(target_right - target_left),
                    word_size(size)
                );
            }
        }
    }

    return r;
}
