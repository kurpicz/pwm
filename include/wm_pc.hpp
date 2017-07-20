/*******************************************************************************
 * include/wm_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once
#ifndef WM_PREFIX_COUNTING
#define WM_PREFIX_COUNTING

#include <cstring>
#include <vector>

#include "common.hpp"
#include "pc.hpp"

template <typename AlphabetType, typename SizeType = uint64_t>
class wm_pc {

public:
    wm_pc() = default;

    wm_pc(const std::vector<AlphabetType>& text,
          const SizeType size,
          const SizeType levels) : _zeros(levels, 0) {

        if(text.size() == 0) { return; }

        _bv = Bvs<SizeType>(size, levels);
        auto ctx = SingleThreaded<SizeType> {
            levels
        };

        SizeType cur_max_char = (1 << levels);


        // While initializing the histogram, we also compute the first level
        SizeType cur_pos = 0;
        for (; cur_pos + 64 <= size; cur_pos += 64) {
            uint64_t word = 0ULL;
            for (SizeType i = 0; i < 64; ++i) {
                ++ctx.hist(levels, text[cur_pos + i]);
                word <<= 1;
                word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
            }
            _bv.vec()[0][cur_pos >> 6] = word;
        }
        if (size & 63ULL) {
            uint64_t word = 0ULL;
            for (SizeType i = 0; i < size - cur_pos; ++i) {
                ++ctx.hist(levels, text[cur_pos + i]);
                word <<= 1;
                word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
            }
            word <<= (64 - (size & 63ULL));
            _bv.vec()[0][size >> 6] = word;
        }

        // The number of 0s at the last level is the number of "even" characters
        ctx.if_zeros([&]() {
            for (SizeType i = 0; i < cur_max_char; i += 2) {
                _zeros[levels - 1] += ctx.hist(levels, i);
            }
        });

        // Now we compute the WM bottom-up, i.e., the last level first
        for (SizeType level = levels - 1; level > 0; --level) {
            const SizeType prefix_shift = (levels - level);
            const SizeType cur_bit_shift = prefix_shift - 1;

            // Update the maximum value of a feasible a bit prefix and update the
            // histogram of the bit prefixes
            cur_max_char >>= 1;
            for (SizeType i = 0; i < cur_max_char; ++i) {
                ctx.hist(level, i)
                    = ctx.hist(level + 1, i << 1) + ctx.hist(level + 1, (i << 1) + 1);
            }

            auto& borders = ctx.borders();

            // Compute the starting positions of characters with respect to their
            // bit prefixes and the bit-reversal permutation
            borders[0] = 0;
            for (SizeType i = 1; i < cur_max_char; ++i) {
                borders[ctx.rho(level, i)] = borders[ctx.rho(level, i - 1)] +
                ctx.hist(level, ctx.rho(level, i - 1));

                ctx.rho(level - 1, i - 1) = (ctx.rho(level, i - 1) >> 1);
            }

            // The number of 0s is the position of the first 1 in the previous level
            ctx.if_zeros([&]() {
               _zeros[level - 1] = borders[1];
            });

            // Now we insert the bits with respect to their bit prefixes
            for (SizeType i = 0; i < size; ++i) {
                const SizeType pos = borders[text[i] >> prefix_shift]++;
                _bv.vec()[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
                << (63ULL - (pos & 63ULL)));
            }
        }
    }

    auto get_bv_and_zeros() const {
        return std::make_pair(_bv.vec(), _zeros);
    }

private:
    Bvs<SizeType> _bv;
    std::vector<SizeType> _zeros;
}; // class wm_pc

#endif // WM_PREFIX_COUNTING

/******************************************************************************/
