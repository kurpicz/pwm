/*******************************************************************************
 * include/util/merge.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/


#pragma once

#include <cassert>
#include <climits>
#include <omp.h>

#include "bit_vectors.hpp"
#include "common.hpp"

template<typename WordType>
void copy_bits(WordType* const dst,
         WordType const* const src,
         uint64_t& dst_off_ref,
         uint64_t& src_off_ref,
         uint64_t const block_size)
{
  if (block_size == 0) return;

  WordType constexpr BITS = (sizeof(WordType) * CHAR_BIT);
  WordType constexpr MOD_MASK = BITS - 1;
  WordType constexpr SHIFT = log2(MOD_MASK);

  uint64_t dst_off = dst_off_ref;
  uint64_t src_off = src_off_ref;
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

      WordType*     ds = dst + (dst_off >> SHIFT);
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

      WordType*     ds = dst + (dst_off >> SHIFT);
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

template<typename ctx_t, typename Rho>
inline auto merge_bit_vectors(uint64_t size,
            uint64_t levels,
            uint64_t shards,
            const std::vector<ctx_t>& src_ctxs,
            const Rho& rho)
{
  assert(shards == src_ctxs.size());

  // Allocate data structures centrally

  struct MergeLevelCtx {
    std::vector<uint64_t> read_offsets;
    uint64_t offset_in_first_word;
    uint64_t first_read_block;
  };

  struct MergeCtx {
    uint64_t offset;
    std::vector<MergeLevelCtx> levels;
  };

  auto ctxs = std::vector<MergeCtx>{shards, {
    0,
    std::vector<MergeLevelCtx> {
      levels, {
        std::vector<uint64_t>(shards),
        0,
        0,
      }
    },
  }};

  // Calculate bit offset per merge shard (thread)
  for (size_t rank = 1; rank < shards; rank++) {
    const size_t omp_rank = rank;
    const size_t omp_size = shards;

    const uint64_t offset = (omp_rank * (word_size(size) / omp_size)) +
      std::min<uint64_t>(omp_rank, word_size(size) % omp_size);

    ctxs[rank - 1].offset = offset * 64ull;
  }
  ctxs[shards - 1].offset = word_size(size) * 64ull;

  //#pragma omp parallel for
  for(size_t level = 0; level < levels; level++) {
    const size_t br_size = 1ull << level;

    size_t write_offset = 0; // bit offset in destination bv
    size_t merge_shard = 0; // index of merge thread

    for(size_t i = 0; i < br_size * shards; i++) {
      const auto bit_i = i / shards;
      const auto read_shard = i % shards;

      auto read_offset =
        [&level, &ctxs](auto merge_shard, auto read_shard) -> uint64_t& {
          return ctxs[merge_shard].levels[level].read_offsets[read_shard];
      };

      auto block_size = src_ctxs[read_shard].hist(level, rho(level, bit_i));

      write_offset += block_size;
      read_offset(merge_shard + 1, read_shard) += block_size;

      // If we passed the current right border, split up the block
      if (write_offset > ctxs[merge_shard].offset) {
        // Take back the last step
        write_offset -= block_size;
        read_offset(merge_shard + 1, read_shard) -= block_size;

        uint64_t offset_in_first_word = 0;
        do {
          // Split up the block like this:
          //     [  left_block_size  |    right_block_size    ]
          //     ^             ^                ^
          // (write_offset)  (ctxs[merge_shard].offset)    (write_offset + block_size)

          auto const left_block_size = ctxs[merge_shard].offset - write_offset;

          write_offset += left_block_size;
          read_offset(merge_shard + 1, read_shard) += left_block_size;

          offset_in_first_word += left_block_size;
          ctxs[merge_shard + 1].levels[level].offset_in_first_word =
            offset_in_first_word;
          ctxs[merge_shard + 1].levels[level].first_read_block = i;

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

  auto r = bit_vectors(size, levels);
  auto& _bv = r.vec();

  #pragma omp parallel
  //for(size_t omp_rank = 0; omp_rank < shards; omp_rank++)
  {
    assert(size_t(omp_get_num_threads()) == shards);
    const size_t omp_rank = omp_get_thread_num();
    const size_t merge_shard = omp_rank;

    auto& ctx = ctxs[merge_shard];

    const auto target_right = std::min(ctx.offset, size);
    const auto target_left = std::min((merge_shard > 0 ?
      ctxs[merge_shard - 1].offset : 0), target_right);

    for (size_t level = 0; level < levels; level++) {
      auto seq_i = ctx.levels[level].first_read_block;

      uint64_t write_offset = target_left;
      size_t init_offset = ctx.levels[level].offset_in_first_word;

      while (write_offset < target_right) {
        const auto i = seq_i / shards;
        const auto read_shard = seq_i % shards;
        seq_i++;

        const auto& h = src_ctxs[read_shard];
        const auto& local_bv = src_ctxs[read_shard].bv().vec()[level];

        auto const block_size = std::min<uint64_t>(
          target_right - write_offset,
          h.hist(level, rho(level, i)) - init_offset
        );
        init_offset = 0; // TODO: remove this by doing a initial pass

        auto& local_cursor = ctx.levels[level].read_offsets[read_shard];

        copy_bits<uint64_t>(
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

/******************************************************************************/
