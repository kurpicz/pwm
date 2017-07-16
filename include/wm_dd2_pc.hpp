/*******************************************************************************
 * include/wm_dd2_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_DOMAIN_DECOMPOSITION2_PREFIX_COUNTING
#define WM_DOMAIN_DECOMPOSITION2_PREFIX_COUNTING

#include <algorithm>
#include <chrono>
#include <cstring>
#include <omp.h>
#include <vector>

#include "common.hpp"
#include "debug.hpp"

template <typename AlphabetType, typename SizeType>
class wm_dd2_pc {

public:
    wm_dd2_pc(const std::vector<AlphabetType>& text,
             const SizeType size,
             const SizeType levels) : _bv(levels), _zeros(levels, 0) {

        if(text.size() == 0) { return; }

        // TODO: constexpr in common.hpp
        auto word_size = [](uint64_t size) {
            return (size + 63ULL) >> 6;
        };

        // TODO: flatten vectors

        const SizeType shards = omp_get_max_threads();

        // This will be used for merging
        auto glob_zeros = std::vector<std::vector<SizeType>>(
            shards, std::vector<SizeType>(levels));
        auto glob_bv = std::vector<std::vector<uint64_t*>>(
            shards
        );
        auto glob_hist = std::vector<std::vector<std::vector<SizeType>>>(
            shards, std::vector<std::vector<SizeType>>(
                levels + 1, std::vector<SizeType>((1 << levels), 0)));

        // This just exists for central allocation
        auto glob_borders = std::vector<std::vector<SizeType>>(
            shards, std::vector<SizeType>((1 << levels), 0));

        #pragma omp parallel
        {
            const size_t omp_rank = omp_get_thread_num();
            const size_t omp_size = omp_get_num_threads();

            const SizeType local_size = (size / omp_size) +
                ((omp_rank < size % omp_size) ? 1 : 0);
            const SizeType offset = (omp_rank * (size / omp_size)) +
                std::min<SizeType>(omp_rank, size % omp_size);

            SizeType cur_max_char = (1 << levels);
            std::vector<SizeType> bit_reverse = BitReverse<SizeType>(levels - 1);

            auto& zeros = glob_zeros[omp_rank];
            auto& borders = glob_borders[omp_rank];
            auto& hist = glob_hist[omp_rank];
            auto& tmp_bv = glob_bv[omp_rank];

            tmp_bv = std::vector<uint64_t*>(levels);

            tmp_bv[0] = new uint64_t[word_size(local_size) * levels];
            // memset is ok (all to 0)
            memset(tmp_bv[0], 0, (word_size(local_size) * sizeof(uint64_t)) * levels);
            for (SizeType level = 1; level < levels; ++level) {
                tmp_bv[level] = tmp_bv[level - 1] + word_size(local_size);
            }

            // While initializing the histogram, we also compute the fist level
            SizeType cur_pos = 0;
            for (; cur_pos + 64 <= local_size; cur_pos += 64) {
                uint64_t word = 0ULL;
                for (SizeType i = offset; i < 64 + offset; ++i) {
                    ++hist[levels][text[cur_pos + i]];
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                tmp_bv[0][cur_pos >> 6] = word;
            }
            if (local_size & 63ULL) {
                uint64_t word = 0ULL;
                for (SizeType i = offset; i < local_size - cur_pos + offset; ++i) {
                    ++hist[levels][text[cur_pos + i]];
                    word <<= 1;
                    word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
                }
                word <<= (64 - (local_size & 63ULL));
                tmp_bv[0][local_size >> 6] = word;
            }

            // The number of 0s at the last level is the number of "even" characters
            for (SizeType i = 0; i < cur_max_char; i += 2) {
                zeros[levels - 1] += hist[levels][i];
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
                    borders[bit_reverse[i]] = borders[bit_reverse[i - 1]] +
                        hist[level][bit_reverse[i - 1]];
                    bit_reverse[i - 1] >>= 1;
                }

                // The number of 0s is the position of the first 1 in the previous level
                zeros[level - 1] = borders[1];

                // Now we insert the bits with respect to their bit prefixes
                for (SizeType i = offset; i < local_size + offset; ++i) {
                    const SizeType pos = borders[text[i] >> prefix_shift]++;
                    tmp_bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
                        << (63ULL - (pos & 63ULL)));
                }
            }

            hist[0].resize(1);
            if (hist.size() > 1) {
                hist[0][0] = hist[1][0] + hist[1][1];
            }
        }

        DCHECK_EQ(shards, glob_bv.size());

        // TODO: Add abstraction for allocating the bitvector (no more bare vector of pointers)

        auto bit_reverse = std::vector<std::vector<SizeType>>(levels);
        bit_reverse[levels - 1] = BitReverse<SizeType>(levels - 1);
        for(size_t level = levels - 1; level > 0; level--) {
            bit_reverse[level - 1] = std::vector<SizeType>(bit_reverse[level].size() / 2);
            for(size_t i = 0; i < bit_reverse[level - 1].size(); i++) {
                bit_reverse[level - 1][i] = bit_reverse[level][i] >> 1;
            }
        }

        auto offsets = std::vector<SizeType>(shards);
        auto local_offsets = std::vector<std::vector<std::vector<SizeType>>>(
            levels, std::vector<std::vector<SizeType>>(
                shards, std::vector<SizeType>(shards + 1)
            )
        );

        auto word_offsets = std::vector<std::vector<SizeType>>(
            levels, std::vector<SizeType>(shards + 1)
        );
        auto block_seq_offsets = std::vector<std::vector<SizeType>>(
            levels, std::vector<SizeType>(shards + 1)
        );

        std::cout << "================================================\n";
        std::cout << "global size: " << size << "\n";
        std::cout << "word size:  " << word_size(size) << "\n";
        std::cout << "================================================\n";
        for (size_t rank = 1; rank < shards; rank++) {
            const size_t omp_rank = rank;
            const size_t omp_size = shards;

            const SizeType offset = (omp_rank * (word_size(size) / omp_size)) +
                std::min<SizeType>(omp_rank, word_size(size) % omp_size);

            offsets[rank - 1] = offset * 64ull;
        }
        offsets[shards - 1] = word_size(size) * 64ull;

        size_t last = 0;
        for (size_t rank = 0; rank < shards; rank++) {
            std::cout << "offset: " << offsets[rank] << ", " << (offsets[rank] - last) << "\n";
            last = offsets[rank];
        }

        std::cout << "================================================\n";

        // TODO: fix
        for(size_t level = 0; level < levels; level++) {
            const auto& br = bit_reverse[level];
            const size_t br_size = 1ull << level;
            size_t j = 0;
            size_t oi = 0;

            auto body = [&](auto j, auto block_size, char sym, size_t x, size_t y) {
                std::cout
                    << sym <<
                    " j: " << std::setw(3)<< (j - block_size)
                    << " += " << std::setw(3)<< block_size
                    << " == " << std::setw(3)<< j
                    << "\n";
                for (size_t a = 0; a < shards; a++) {
                    std::cout <<
                        "  shard: " << a <<
                        ", offset: ";
                    for (size_t b = 0; b < shards; b++) {
                        char mark = ' ';
                        if (a == y && b == x) mark = '+';
                        std::cout << mark << std::setw(3) << local_offsets[level][a][b] << ", ";
                    }
                    char mark = ' ';
                    if (a == y && shards == x) mark = '+';
                    std::cout << "| " << mark << std::setw(3) << local_offsets[level][a][shards];
                    std::cout << "\n";
                }
            };

            for(size_t i = 0; i < br_size * shards; i++) {
                const auto bit_i = i / shards;
                const auto shard = i % shards;

                const auto h = glob_hist[shard][level];
                const auto block_size = h[br[bit_i]];

                j += block_size;
                local_offsets[level][shard][oi + 1] += block_size;

                if (j <= offsets[oi]) {
                    body(j, block_size, 'm', oi + 1, shard);
                } else {
                    size_t right = j - offsets[oi];
                    size_t left = block_size - right;

                    j -= block_size;
                    local_offsets[level][shard][oi + 1] -= block_size;

                    j += left;
                    local_offsets[level][shard][oi + 1] += left;
                    body(j, left, 'l', oi + 1, shard);

                    // TODO: rename word_offset to something like
                    // "offset in block"
                    word_offsets[level][oi + 1] = left;
                    block_seq_offsets[level][oi + 1] = i;

                    std::cout << "------------------------------------------------\n";
                    if (oi + 2 < shards + 1) {
                        for(size_t s = 0; s < shards; s++) {
                            local_offsets[level][s][oi + 2] = local_offsets[level][s][oi + 1];
                        }
                    }
                    oi++;

                    // TODO: dcheck about block_seq_offsets, word_offsets

                    j += right;
                    local_offsets[level][shard][oi + 1] += right;
                    body(j, right, 'r', oi + 1, shard);
                }

                //body(j, 0, 'f');
            }
            //std::cout << "j: " << j << "\n";
            std::cout << "================================================\n";
        }

        // TODO: print out all data gathered before

        for (size_t level = 0; level < levels; level++) {
            for(size_t merge_shard = 0; merge_shard < shards; merge_shard++) {
                std::cout << "merge shard " << merge_shard << ": ";

                std::cout << "offsets = [";
                for(size_t read_shard = 0; read_shard < shards; read_shard++) {
                    std::cout << std::setw(3) << local_offsets[level][read_shard][merge_shard] << ", ";
                }
                std::cout << "]";

                std::cout << ", block_seq_offset: " << block_seq_offsets[level][merge_shard];
                std::cout << "(" << block_seq_offsets[level][merge_shard] / shards;
                std::cout << ", " << block_seq_offsets[level][merge_shard] % shards;
                std::cout << "), block_offset: " << word_offsets[level][merge_shard];

                std::cout << "\n";
            }

            std::cout << "\n";
        }

        std::cout << "================================================\n";

        // TODO: factor out
        for(size_t level = 0; level < levels; level++) {
            _bv[level] = new uint64_t[(size + 63ULL) >> 6];
            memset(_bv[level], 0, word_size(size) * sizeof(uint64_t));
        }

        //#pragma omp parallel
        //{
        //const size_t omp_rank = omp_get_thread_num();
        //const size_t omp_size = omp_get_num_threads();

        //#pragma omp parallel for
        for(size_t merge_shard = 0; merge_shard < shards; merge_shard++) {
            const auto target_right = std::min(offsets[merge_shard], size);
            const auto target_left = std::min((merge_shard > 0 ? offsets[merge_shard - 1] : 0), target_right);
            const auto target_size = target_right - target_left;

            auto cursors = std::vector<std::vector<size_t>>(
                levels, std::vector<size_t>(shards)
            );

            for (size_t level = 0; level < levels; level++) {
                for(size_t read_shard = 0; read_shard < shards; read_shard++) {
                    cursors[level][read_shard] = local_offsets[level][read_shard][merge_shard];
                }
            }

            for (size_t level = 0; level < levels; level++) {
                const size_t br_size = 1ull << level;
                const auto& br = bit_reverse[level];

                auto seq_i = block_seq_offsets[level][merge_shard];

                {
                    const auto i = seq_i / shards;
                    const auto shard = seq_i % shards;

                    std::cout << "read starting @ block "
                        << i << "," << shard
                        <<  " so many bits: " << target_size
                        << " (" << target_left << " - " << target_right << ")"
                        << " with initial offset " <<  word_offsets[level][merge_shard];
                    std::cout << ", cursors = [";
                    for(size_t read_shard = 0; read_shard < shards; read_shard++) {
                        std::cout << std::setw(3) << cursors[level][read_shard] << ", ";
                    }
                    std::cout << "]\n";
                }


                size_t j = target_left;
                size_t init_offset = word_offsets[level][merge_shard];

                while (true) {
                    const auto i = seq_i / shards;
                    const auto shard = seq_i % shards;

                    const auto& h = glob_hist[shard][level];
                    const auto& local_bv = glob_bv[shard][level];

                    auto block_size = h[br[i]] - init_offset;
                    init_offset = 0; // TODO: remove this by doing a initial pass
                    auto& local_cursor = cursors[level][shard];

                    // TODO: copy over whole block
                    while(block_size != 0) {
                        block_size--;
                        const auto src_pos = local_cursor++;
                        const auto pos = j++;
                        const bool bit = bit_at(local_bv, src_pos);

                        _bv[level][pos >> 6] |= (uint64_t(bit) << (63ULL - (pos & 63ULL)));

                        if (j >= target_right) {
                            break;
                        }
                    }

                    if (j >= target_right) {
                        break;
                    }

                    seq_i++;
                }



            }
            std::cout << "\n";
        }

        auto cursors = std::vector<std::vector<size_t>>(
            levels, std::vector<size_t>(shards)
        );

        //#pragma omp parallel for
        for(size_t level = 0; level < levels; level++) {
            const auto& br = bit_reverse[level];
            size_t j = 0;
            const size_t br_size = 1ull << level;

            for(size_t seq_i = 0; seq_i < br_size * shards; seq_i++) {
                const auto i = seq_i / shards;
                const auto shard = seq_i % shards;

                const auto& h = glob_hist[shard][level];
                const auto& local_bv = glob_bv[shard][level];

                auto block_size = h[br[i]];
                auto& local_cursor = cursors[level][shard];

                // TODO: copy over whole block
                while(block_size != 0) {
                    block_size--;
                    const auto src_pos = local_cursor++;
                    const auto pos = j++;
                    const bool bit = bit_at(local_bv, src_pos);

                    //_bv[level][pos >> 6] |= (uint64_t(bit) << (63ULL - (pos & 63ULL)));
                }
            }
            for(size_t shard = 0; shard < glob_bv.size(); shard++) {
                _zeros[level] += glob_zeros[shard][level];
            }
        }

        std::cout << "single merged:\n";
        print_bv(_bv, size);
    }

    auto get_bv_and_zeros() const {
        return std::make_pair(_bv, _zeros);
    }

private:
    std::vector<uint64_t*> _bv;
    std::vector<SizeType> _zeros;
}; // class wm_dd2_pc

#endif // WM_DOMAIN_DECOMPOSITION_PREFIX_COUNTING

/******************************************************************************/
