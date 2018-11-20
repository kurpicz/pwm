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

#include "arrays/bit_vectors.hpp"
#include "util/common.hpp"
#include "util/macros.hpp"

template <typename WordType, bool zero_mem = false>
void copy_bits(WordType* const dst,
               WordType const* const src,
               uint64_t& dst_off_ref,
               uint64_t& src_off_ref,
               uint64_t const block_size,
               WordType const** opt_zero_mem_marker = nullptr) {
  if (block_size == 0) {
    return;
  }

  /*
  NB: Code skeleton for debugging writes to memory.
  #include <bitset>
  #include <iomanip>
  #include <sstream>
  auto debug_ptr = [dst](WordType* const ptr, auto extra_fmt) {
    int64_t diff = ptr - dst;
    std::stringstream ss;
    ss << "write [" << dst << "][" << std::setw(6) << diff
       << "]: " << std::bitset<64>(*ptr);
    extra_fmt(ss);
    ss << "\n";
    // Using a stringstream to get a atomic write to cout
    std::cout << ss.str();
  };
  */

  auto zero_word_once = [opt_zero_mem_marker](WordType* word) {
    if constexpr (zero_mem) {
      if (opt_zero_mem_marker != nullptr) {
        WordType const*& zero_mem_marker = *opt_zero_mem_marker;
        if (zero_mem_marker == nullptr) {
          zero_mem_marker = word;
        }
        if (zero_mem_marker == word) {
          *word = 0;
          zero_mem_marker++;
        }
        DCHECK(zero_mem_marker == word + 1);
      }
    } else {
      (void) opt_zero_mem_marker;
      (void) word;
    }
  };
  auto skip_zero_words = [opt_zero_mem_marker](WordType const* start,
                                               size_t count) {
    if constexpr (zero_mem) {
      if (opt_zero_mem_marker != nullptr) {
        WordType const*& zero_mem_marker = *opt_zero_mem_marker;
        if (zero_mem_marker == nullptr) {
          zero_mem_marker = start;
        }
        zero_mem_marker += count;
      }
    } else {
      (void) opt_zero_mem_marker;
      (void) start;
      (void) count;
    }
  };

  WordType constexpr BITS = (sizeof(WordType) * CHAR_BIT);
  WordType constexpr MOD_MASK = BITS - 1;
  WordType constexpr SHIFT = pwmlog2(MOD_MASK);

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

      if ((dst_off & MOD_MASK) != 0 && dst_off != dst_off_end) {
        zero_word_once(&word);
      }
      while ((dst_off & MOD_MASK) != 0 && dst_off != dst_off_end) {
        bool const bit = bit_at<WordType>(src, src_off++);

        word |= (WordType(bit) << (MOD_MASK - (dst_off++ & MOD_MASK)));
      }
    }

    // Copy the the bulk in-between word-wise
    {
      auto const words = (dst_off_end - dst_off) >> SHIFT;

      WordType* ds = dst + (dst_off >> SHIFT);
      WordType const* sr = src + (src_off >> SHIFT);

      WordType const* const ds_end = ds + words;

      skip_zero_words(ds, words);
      while (ds != ds_end) {
        *ds++ = *sr++;
      }

      dst_off += words * BITS;
      src_off += words * BITS;
    }

    // Copy unaligned trailing bits
    {
      auto& word = dst[dst_off >> SHIFT];

      if (dst_off != dst_off_end) {
        zero_word_once(&word);
      }
      while (dst_off != dst_off_end) {
        bool const bit = bit_at<WordType>(src, src_off++);

        word |= (WordType(bit) << (MOD_MASK - (dst_off++ & MOD_MASK)));
      }
    }
  } else {
    // Copy unaligned leading bits
    {
      auto& word = dst[dst_off >> SHIFT];

      if ((dst_off & MOD_MASK) != 0 && dst_off != dst_off_end) {
        zero_word_once(&word);
      }
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

      WordType* ds = dst + (dst_off >> SHIFT);
      WordType const* sr = src + (src_off >> SHIFT);
      WordType const* const ds_end = ds + words;

      skip_zero_words(ds, words);
      while (ds != ds_end) {
        *ds++ = (*sr << src_shift_a) | (*(sr + 1) >> src_shift_b);
        sr++;
      }

      dst_off += words * BITS;
      src_off += words * BITS;
    }

    // Copy unaligned trailing bits
    {
      auto& word = dst[dst_off >> SHIFT];

      if (dst_off != dst_off_end) {
        zero_word_once(&word);
      }
      while (dst_off != dst_off_end) {
        bool const bit = bit_at<WordType>(src, src_off++);

        word |= (WordType(bit) << (MOD_MASK - (dst_off++ & MOD_MASK)));
      }
    }
  }

  dst_off_ref += block_size;
  src_off_ref += block_size;
}

template <typename ContextType, typename Rho>
inline auto merge_bit_vectors(uint64_t size,
                              uint64_t levels,
                              uint64_t shards,
                              const std::vector<ContextType>& src_ctxs,
                              const Rho& rho) {
  assert(shards == src_ctxs.size());

  // Allocate data structures centrally

  struct MergeLevelCtx {
    std::vector<uint64_t> read_offsets;
    uint64_t initial_offset;
    uint64_t first_read_block;
  };

  struct MergeCtx {
    uint64_t end_offset;
    std::vector<MergeLevelCtx> levels;
  };

  auto ctxs = std::vector<MergeCtx>{
      shards,
      {
          0,
          std::vector<MergeLevelCtx>{levels,
                                     {
                                         std::vector<uint64_t>(shards),
                                         0,
                                         0,
                                     }},
      }};

  /*
    Visualization of ctxs for levels = 2 and shards = 2:

    ┌──────────────┬───────────────────────┐
    │ctxs: MergeCtx│ x shards              │
    │ ┌────────────┴───────────────────────┤
    │ │end_offset: uint64_t                │
    │ ├─────────────────────┬──────────────┤
    │ │levels: MergeLevelCtx│ x levels     │
    │ │ ┌───────────────────┴──────────────┤
    │ │ │initial_offset: uint64_t          │
    │ │ ├──────────────────────────────────┤
    │ │ │first_read_block: uint64_t        │
    │ │ ├──────────────────────┬───────────┤
    │ │ │read_offsets: uint64_t│ x shards  │
    │ │ │ ┌────────────────────┴───────────┤
    │ │ │ │                                │
    │ │ │ ╞════════════════════════════════╡
    │ │ │ │                                │
    │ │ ╞═╧════════════════════════════════╡
    │ │ │initial_offset: uint64_t          │
    │ │ ├──────────────────────────────────┤
    │ │ │first_read_block: uint64_t        │
    │ │ ├──────────────────────┬───────────┤
    │ │ │read_offsets: uint64_t│ x shards  │
    │ │ │ ┌────────────────────┴───────────┤
    │ │ │ │                                │
    │ │ │ ╞════════════════════════════════╡
    │ │ │ │                                │
    │ ╞═╧═╧════════════════════════════════╡
    │ │end_offset: uint64_t                │
    │ ├─────────────────────┬──────────────┤
    │ │levels: MergeLevelCtx│ x levels     │
    │ │ ┌───────────────────┴──────────────┤
    │ │ │initial_offset: uint64_t          │
    │ │ ├──────────────────────────────────┤
    │ │ │first_read_block: uint64_t        │
    │ │ ├──────────────────────┬───────────┤
    │ │ │read_offsets: uint64_t│ x shards  │
    │ │ │ ┌────────────────────┴───────────┤
    │ │ │ │                                │
    │ │ │ ╞════════════════════════════════╡
    │ │ │ │                                │
    │ │ ╞═╧════════════════════════════════╡
    │ │ │initial_offset: uint64_t          │
    │ │ ├──────────────────────────────────┤
    │ │ │first_read_block: uint64_t        │
    │ │ ├──────────────────────┬───────────┤
    │ │ │read_offsets: uint64_t│ x shards  │
    │ │ │ ┌────────────────────┴───────────┤
    │ │ │ │                                │
    │ │ │ ╞════════════════════════════════╡
    │ │ │ │                                │
    └─┴─┴─┴────────────────────────────────┘

  */

  // Calculate end bit offset per merge shard (thread).
  // It is an multiple of 64, to ensure no interference
  // in writing bits to the left or right of it in parallel
  for (size_t rank = 1; rank < shards; rank++) {
    const size_t omp_rank = rank;
    const size_t omp_size = shards;

    const uint64_t offset =
        (omp_rank * (word_size(size) / omp_size)) +
        std::min<uint64_t>(omp_rank, word_size(size) % omp_size);

    ctxs[rank - 1].end_offset = offset * 64ull;
  }
  ctxs[shards - 1].end_offset = word_size(size) * 64ull;

  //#pragma omp parallel for
  for (size_t level = 0; level < levels; level++) {
    // number of tree nodes on the current level
    const size_t blocks = 1ull << level;

    size_t write_offset = 0; // bit offset in destination bv
    size_t merge_shard = 0;  // index of merge thread

    // iterate over all blocks of all shards on the current level
    //
    // this gradually assigns the bits of the
    // blocks to each of the aviable merge shards of the current level:
    //
    // |   read_shard  |   read_shard  |
    // +---------------+---------------+
    // |     block     |     block     |
    // | block | block | block | block | <- current level
    // +----------+----------+---------+
    // | m. shard | m. shard |m. shard |
    //
    for (size_t i = 0; i < blocks * shards; i++) {
      // returns merge level context of merge_shard+1
      auto nxt_lctx = [&level, &ctxs](auto merge_shard) -> MergeLevelCtx& {
        return ctxs[merge_shard + 1].levels[level];
      };

      // which block (node on current level of tree)
      const auto block = i / shards;

      // which shard to read from
      const auto read_shard = i % shards;

      // map logical block index to
      // actual block index (rho depends on WT vs WM)
      const auto permuted_block = rho(level, block);

      // block size == number of entries in the block on this level
      auto block_size = src_ctxs[read_shard].hist_at_level(level)[permuted_block];

      // advance global write offset by the number of bits assigned for
      // this block
      write_offset += block_size;

      // advance local read offset of next merge shard
      // for the curent read shard
      //
      // this way, all merge shard start reading
      // at different offsets from the read shards
      nxt_lctx(merge_shard).read_offsets[read_shard] += block_size;

      // If we passed the current right border, split up the block
      if (write_offset > ctxs[merge_shard].end_offset) {
        // Take back the last step
        write_offset -= block_size;
        nxt_lctx(merge_shard).read_offsets[read_shard] -= block_size;

        uint64_t initial_offset = 0;
        do {
          // Split up the block like this:
          //     [  left_block_size  |    right_block_size    ]
          //     ^                   ^                        ^
          // (write_offset)          |                        |
          //           (ctxs[merge_shard].end_offset)         |
          //                                     (write_offset + block_size)
          //
          // this loop iterates multiple times, to ensure right_block_size
          // did not also overlap the next end_offset

          auto const left_block_size =
              ctxs[merge_shard].end_offset - write_offset;

          // advance global and local read offsets to end exactly at
          // end_offset
          write_offset += left_block_size;
          nxt_lctx(merge_shard).read_offsets[read_shard] += left_block_size;

          initial_offset += left_block_size;

          // from which offset in which block to start reading
          nxt_lctx(merge_shard).initial_offset = initial_offset;
          nxt_lctx(merge_shard).first_read_block = i;

          // if there is a merge_shard to the right
          // of the next merge shard, intialize
          // its local read offsets with those of the next shard
          if (merge_shard + 2 < shards) {
            for (size_t s = 0; s < shards; s++) {
              nxt_lctx(merge_shard + 1).read_offsets[s] =
                  nxt_lctx(merge_shard).read_offsets[s];
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
        } while ((write_offset + block_size) > ctxs[merge_shard].end_offset);

        // Process remainder of block
        write_offset += block_size;
        nxt_lctx(merge_shard).read_offsets[read_shard] += block_size;

        assert(write_offset <= ctxs[merge_shard].end_offset);
      }
    }
  triple_loop_exit:; // we are done
  }

  auto bv = bit_vectors<false>(levels, size);

  #pragma omp parallel for
  for (size_t merge_shard = 0; merge_shard < shards; merge_shard++) {
    assert(size_t(omp_get_num_threads()) == shards);

    auto& ctx = ctxs[merge_shard];

    const auto target_right = std::min(ctx.end_offset, size);
    const auto target_left = std::min(
        (merge_shard > 0 ? ctxs[merge_shard - 1].end_offset : 0), target_right);

    for (size_t level = 0; level < levels; level++) {
      auto i = ctx.levels[level].first_read_block;
      uint64_t write_offset = target_left;
      uint64_t const* zero_marker = nullptr;

      auto copy_next_block = [&](size_t const initial_offset) {
        const auto block = i / shards;
        const auto read_shard = i % shards;
        i++;

        const auto& local_bv = src_ctxs[read_shard].bv()[level];
        auto&& hist = src_ctxs[read_shard].hist_at_level(level);
        uint64_t block_size = hist[rho(level, block)] - initial_offset;
        uint64_t distance_to_end = target_right - write_offset;

        uint64_t copy_size;
        if (PWM_LIKELY(block_size <= distance_to_end)) {
          copy_size = block_size;
        } else {
          copy_size = distance_to_end;
        }

        auto& local_cursor = ctx.levels[level].read_offsets[read_shard];
        copy_bits<uint64_t, true>(bv[level].data(), local_bv.data(),
                                  write_offset, local_cursor, copy_size,
                                  &zero_marker);
      };

      if (write_offset < target_right) {
        // the first block might start somewhere in the middle
        copy_next_block(ctx.levels[level].initial_offset);
      }
      while (write_offset < target_right) {
        // other blocks start at 0
        copy_next_block(0);
      }
    }
  }

  return bv;
}

/******************************************************************************/
