#pragma once

#include "common.hpp"
#include "context.hpp"

template<
    typename Text,
    typename SizeType,
    typename ctx_t
>
void pc(Text const& text,
        SizeType const size,
        SizeType const levels,
        ctx_t& ctx)
{
    SizeType cur_max_char = (1 << levels);

    auto& zeros = ctx.zeros();
    auto& borders = ctx.borders();
    auto& bv = ctx.bv().vec();

    // While initializing the histogram, we also compute the first level
    SizeType cur_pos = 0;
    for (; cur_pos + 64 <= size; cur_pos += 64) {
        uint64_t word = 0ULL;
        for (SizeType i = 0; i < 64; ++i) {
            ++ctx.hist(levels, text[cur_pos + i]);
            word <<= 1;
            word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        bv[0][cur_pos >> 6] = word;
    }
    if (size & 63ULL) {
        uint64_t word = 0ULL;
        for (SizeType i = 0; i < size - cur_pos; ++i) {
            ++ctx.hist(levels, text[cur_pos + i]);
            word <<= 1;
            word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        bv[0][size >> 6] = word;
    }

    // The number of 0s at the last level is the number of "even" characters
    if (ctx_t::compute_zeros) {
        for (SizeType i = 0; i < cur_max_char; i += 2) {
            zeros[levels - 1] += ctx.hist(levels, i);
        }
    }

    // Now we compute the WM bottom-up, i.e., the last level first
    for (SizeType level = levels - 1; level > 0; --level) {
        const SizeType prefix_shift = (levels - level);
        const SizeType cur_bit_shift = prefix_shift - 1;

        // Update the maximum value of a feasible a bit prefix and update the
        // histogram of the bit prefixes
        cur_max_char >>= 1;
        for (SizeType i = 0; i < cur_max_char; ++i) {
            ctx.hist(level, i)
                = ctx.hist(level + 1, i << 1)
                + ctx.hist(level + 1, (i << 1) + 1);
        }

        // Compute the starting positions of characters with respect to their
        // bit prefixes and the bit-reversal permutation
        borders[0] = 0;
        for (SizeType i = 1; i < cur_max_char; ++i) {
            auto const prev_rho = ctx.rho(level, i - 1);

            borders[ctx.rho(level, i)] = borders[prev_rho] + ctx.hist(level, prev_rho);

            if (ctx_t::compute_rho)  {
                ctx.set_rho(level - 1, i - 1, prev_rho >> 1);
            }
        }

        // The number of 0s is the position of the first 1 in the previous level
        if (ctx_t::compute_zeros) {
            zeros[level - 1] = borders[1];
        }

        // Now we insert the bits with respect to their bit prefixes
        for (SizeType i = 0; i < size; ++i) {
            const SizeType pos = borders[text[i] >> prefix_shift]++;
            bv[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
            << (63ULL - (pos & 63ULL)));
        }
    }

    if (levels > 1) { // TODO check condition
        ctx.hist(0, 0) = ctx.hist(1, 0) + ctx.hist(1, 1);
    }
}
