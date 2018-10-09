/*******************************************************************************
 * include/huffman/huff_merge.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cassert>
#include <climits>
#include <omp.h>

#include "construction/merge.hpp"
#include "huffman/huff_bit_vectors.hpp"
#include "util/common.hpp"
#include "util/macros.hpp"

template <typename ContextType, typename Rho>
inline auto huff_merge_bit_vectors(std::vector<uint64_t> const& level_sizes,
                                   uint64_t const shards,
                                   const std::vector<ContextType>& src_ctxs,
                                   const Rho& rho) {
  assert(shards == src_ctxs.size());
  assert(shards > 1);

  uint64_t const levels = level_sizes.size();

  // Allocate data structures centrally

  struct MergeLevelCtx {
    std::vector<uint64_t> read_offsets;
    uint64_t initial_read_offset;
    uint64_t first_read_block;
    uint64_t write_end_offset;
  };

  struct MergeCtx {
    std::vector<MergeLevelCtx> levels;
  };

  auto ctxs = std::vector<MergeCtx>{
      shards,
      {
          std::vector<MergeLevelCtx>{levels,
                                     {
                                         std::vector<uint64_t>(shards),
                                         0,
                                         0,
                                         0,
                                     }},
      }};

  /*
    Visualization of ctxs for levels = 2 and shards = 2:
    ┌──────────────┬───────────────────────┐
    │ctxs: MergeCtx│ x shards              │
    │ ┌────────────┴────────┬──────────────┤
    │ │levels: MergeLevelCtx│ x levels     │
    │ │ ┌───────────────────┴──────────────┤
    │ │ │write_end_offset: uint64_t        │
    │ │ ├──────────────────────────────────┤
    │ │ │initial_read_offset: uint64_t     │
    │ │ ├──────────────────────────────────┤
    │ │ │first_read_block: uint64_t        │
    │ │ ├──────────────────────┬───────────┤
    │ │ │read_offsets: uint64_t│ x shards  │
    │ │ │ ┌────────────────────┴───────────┤
    │ │ │ │                                │
    │ │ │ ╞════════════════════════════════╡
    │ │ │ │                                │
    │ │ ╞═╧════════════════════════════════╡
    │ │ │initial_read_offset: uint64_t     │
    │ │ ├──────────────────────────────────┤
    │ │ │first_read_block: uint64_t        │
    │ │ ├──────────────────────┬───────────┤
    │ │ │read_offsets: uint64_t│ x shards  │
    │ │ │ ┌────────────────────┴───────────┤
    │ │ │ │                                │
    │ │ │ ╞════════════════════════════════╡
    │ │ │ │                                │
    │ ╞═╧═╧════════════════════════════════╡
    │ │levels: MergeLevelCtx│ x levels     │
    │ │ ┌───────────────────┴──────────────┤
    │ │ │write_end_offset: uint64_t        │
    │ │ ├──────────────────────────────────┤
    │ │ │initial_read_offset: uint64_t     │
    │ │ ├──────────────────────────────────┤
    │ │ │first_read_block: uint64_t        │
    │ │ ├──────────────────────┬───────────┤
    │ │ │read_offsets: uint64_t│ x shards  │
    │ │ │ ┌────────────────────┴───────────┤
    │ │ │ │                                │
    │ │ │ ╞════════════════════════════════╡
    │ │ │ │                                │
    │ │ ╞═╧════════════════════════════════╡
    │ │ │initial_read_offset: uint64_t     │
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
  for (size_t level = 0; level < levels; level++) {
    const size_t size = level_sizes[level];
    for (size_t shard = 1; shard < shards; shard++) {
      const uint64_t offset =
          (shard * (word_size(size) / shards)) +
          std::min<uint64_t>(shard, word_size(size) % shards);

      ctxs[shard - 1].levels[level].write_end_offset =
          std::min<uint64_t>(offset * 64ull, size);
    }
    ctxs[shards - 1].levels[level].write_end_offset =
        std::min<uint64_t>(word_size(size) * 64ull, size);
  }

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
      // std::cout << "[offsets]   i: " << i << "\n";

      // returns merge level context of merge_shard
      auto lctx = [&level, &ctxs](auto merge_shard) -> MergeLevelCtx& {
        return ctxs[merge_shard].levels[level];
      };

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
      auto block_size = src_ctxs[read_shard].hist(level, permuted_block);

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
      if (write_offset > lctx(merge_shard).write_end_offset) {
        // Take back the last step
        write_offset -= block_size;
        nxt_lctx(merge_shard).read_offsets[read_shard] -= block_size;

        uint64_t initial_read_offset = 0;
        do {
          // Split up the block like this:
          //     [  left_block_size  |    right_block_size    ]
          //     ^                   ^                        ^
          // (write_offset)          |                        |
          //           (ctxs[merge_shard].write_end_offset)   |
          //                              (write_offset + block_size)
          //
          // this loop iterates multiple times, to ensure right_block_size
          // did not also overlap the next write_end_offset

          auto const left_block_size =
              lctx(merge_shard).write_end_offset - write_offset;

          // advance global and local read offsets to end exactly at
          // write_end_offset
          write_offset += left_block_size;
          nxt_lctx(merge_shard).read_offsets[read_shard] += left_block_size;

          initial_read_offset += left_block_size;

          // from which offset in which block to start reading
          nxt_lctx(merge_shard).initial_read_offset = initial_read_offset;
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
        } while ((write_offset + block_size) >
                 lctx(merge_shard).write_end_offset);

        // Process remainder of block
        write_offset += block_size;
        nxt_lctx(merge_shard).read_offsets[read_shard] += block_size;

        assert(write_offset <= lctx(merge_shard).write_end_offset);
      }
    }
  triple_loop_exit:; // we are done
  }

  // TODO: remove redundant argument
  auto r = huff_bit_vectors(levels, level_sizes);
  auto& _bv = r;

#pragma omp parallel for
  for (size_t merge_shard = 0; merge_shard < shards; merge_shard++) {
    for (size_t level = 0; level < levels; level++) {
      auto& lctx = ctxs[merge_shard].levels[level];

      const auto target_right =
          std::min(lctx.write_end_offset, level_sizes[level]);
      const auto target_left =
          std::min((merge_shard > 0
                        ? ctxs[merge_shard - 1].levels[level].write_end_offset
                        : 0),
                   target_right);

      uint64_t i = lctx.first_read_block;
      uint64_t write_offset = target_left;

      auto copy_next_block = [&](size_t const initial_read_offset) {
        const auto block = i / shards;
        const auto read_shard = i % shards;
        i++;

        const auto& local_bv = src_ctxs[read_shard].bv()[level];
        const auto& h = src_ctxs[read_shard];
        auto const permuted_block = rho(level, block);
        uint64_t block_size =
            h.hist(level, permuted_block) - initial_read_offset;
        uint64_t distance_to_end = target_right - write_offset;

        uint64_t copy_size;
        if (PWM_LIKELY(block_size <= distance_to_end)) {
          copy_size = block_size;
        } else {
          copy_size = distance_to_end;
        }

        auto& local_cursor = lctx.read_offsets[read_shard];

        copy_bits<uint64_t>(_bv[level].data(), local_bv.data(), write_offset,
                            local_cursor, copy_size);
      };

      if (write_offset < target_right) {
        // the first block might start somewhere in the middle
        copy_next_block(lctx.initial_read_offset);
      }
      while (write_offset < target_right) {
        // other blocks start at 0
        copy_next_block(0);
      }
      assert(write_offset == target_right);
    }
  }

  return r;
}

/******************************************************************************/
