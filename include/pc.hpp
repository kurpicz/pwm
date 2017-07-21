#pragma once

#include "common.hpp"

// TODO: WM/WT abstract that selects zeros and rho

// Overwrite information for each level
template<typename SizeType, bool is_matrix>
struct LevelSinglePass {
    std::vector<SizeType> m_hist;
    std::vector<SizeType> m_borders;
    std::vector<SizeType> m_zeros;
    Bvs<SizeType> m_bv;
    std::vector<SizeType> m_bit_reverse;

    LevelSinglePass(SizeType const size, SizeType const levels) {
        auto cur_max_char = (1ull << levels);

        m_hist.reserve(cur_max_char);
        m_hist.resize(cur_max_char, 0);

        if (is_matrix) {
            m_bit_reverse = BitReverse<SizeType>(levels - 1);
        }

        m_borders.reserve(cur_max_char);
        m_borders.resize(cur_max_char, 0);

        m_zeros.reserve(levels);
        m_zeros.resize(levels, 0);

        m_bv = Bvs<SizeType>(size, levels);
    }

    SizeType const hist_size(SizeType const level) {
        return 1ull << level;
    }

    SizeType& hist(SizeType const level, SizeType const i) {
        return m_hist[i];
    }

    SizeType rho(size_t level, size_t i) {
        if (is_matrix) {
            return m_bit_reverse[i];
        }else {
            return i;
        }
    }

    void set_rho(size_t level, size_t i, SizeType val) {
        if (is_matrix) {
            m_bit_reverse[i] = val;
        }
    }

    std::vector<SizeType>& borders() {
        return m_borders;
    }

    bool calc_zeros() {
        return is_matrix;
    }

    bool calc_rho() {
        return is_matrix;
    }

    std::vector<SizeType>& zeros() {
        return m_zeros;
    }

    Bvs<SizeType>& bv() {
        return m_bv;
    }
};

/// Keep calculated information for individual levels around
template<typename SizeType, bool is_matrix>
struct KeepLevel {
    std::vector<std::vector<SizeType>> m_hist;
    decltype(rho_bit_reverse(0)) const* m_rho;
    std::vector<SizeType> m_borders;
    std::vector<SizeType> m_zeros;
    Bvs<SizeType> m_bv;

    KeepLevel(SizeType const size, SizeType const levels, decltype(rho_bit_reverse(0)) const& rho) {
        auto cur_max_char = (1ull << levels);

        m_hist.reserve(levels + 1);
        m_hist.resize(levels + 1);

        for(size_t level = 0; level < (levels + 1); level++) {
            m_hist[level].reserve(hist_size(level));
            m_hist[level].resize(hist_size(level));
        }

        m_rho = &rho;

        m_borders.reserve(cur_max_char);
        m_borders.resize(cur_max_char, 0);

        m_zeros.reserve(levels);
        m_zeros.resize(levels, 0);

        m_bv = Bvs<SizeType>(size, levels);
    }

    SizeType const hist_size(SizeType const level) {
        return 1ull << level;
    }

    SizeType& hist(SizeType const level, SizeType const i) {
        return m_hist[level][i];
    }

    SizeType& rho(size_t level, size_t i) {
        return (*m_rho)(level, i);
    }

    void set_rho(size_t level, size_t i, SizeType val) {
        // m_rho is already calculated
    }

    std::vector<SizeType>& borders() {
        return m_borders;
    }

    bool calc_zeros() {
        return is_matrix;
    }

    bool calc_rho() {
        return false;
    }

    std::vector<SizeType>& zeros() {
        return m_zeros;
    }

    Bvs<SizeType>& bv() {
        return m_bv;
    }
};


template<
    typename Text,
    typename SizeType,
    typename Context
>
void pc(Text const& text,
        SizeType const size,
        SizeType const levels,
        Context& ctx)
{
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
        ctx.bv().vec()[0][cur_pos >> 6] = word;
    }
    if (size & 63ULL) {
        uint64_t word = 0ULL;
        for (SizeType i = 0; i < size - cur_pos; ++i) {
            ++ctx.hist(levels, text[cur_pos + i]);
            word <<= 1;
            word |= ((text[cur_pos + i] >> (levels - 1)) & 1ULL);
        }
        word <<= (64 - (size & 63ULL));
        ctx.bv().vec()[0][size >> 6] = word;
    }

    // The number of 0s at the last level is the number of "even" characters
    if (ctx.calc_zeros()) {
        for (SizeType i = 0; i < cur_max_char; i += 2) {
            ctx.zeros()[levels - 1] += ctx.hist(levels, i);
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
                = ctx.hist(level + 1, i << 1) + ctx.hist(level + 1, (i << 1) + 1);
        }

        auto& borders = ctx.borders();

        // Compute the starting positions of characters with respect to their
        // bit prefixes and the bit-reversal permutation
        borders[0] = 0;
        for (SizeType i = 1; i < cur_max_char; ++i) {
            auto const prev_rho = ctx.rho(level, i - 1);

            borders[ctx.rho(level, i)] = borders[prev_rho] + ctx.hist(level, prev_rho);

            if (ctx.calc_rho())  {
                ctx.set_rho(level - 1, i - 1, prev_rho >> 1);
            }
        }

        // The number of 0s is the position of the first 1 in the previous level
        if (ctx.calc_zeros()) {
            ctx.zeros()[level - 1] = borders[1];
        }

        // Now we insert the bits with respect to their bit prefixes
        for (SizeType i = 0; i < size; ++i) {
            const SizeType pos = borders[text[i] >> prefix_shift]++;
            ctx.bv().vec()[level][pos >> 6] |= (((text[i] >> cur_bit_shift) & 1ULL)
            << (63ULL - (pos & 63ULL)));
        }
    }
}
