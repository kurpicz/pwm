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

    static bool constexpr compute_zeros = is_matrix;
    static bool constexpr compute_rho = is_matrix;

    std::vector<SizeType>& zeros() {
        return m_zeros;
    }

    Bvs<SizeType>& bv() {
        return m_bv;
    }
};

/// Keep calculated information for individual levels around
template<typename SizeType, bool is_matrix, typename rho_t>
struct KeepLevel {
    std::vector<std::vector<SizeType>> m_hist;
    rho_t const* m_rho = nullptr;
    std::vector<SizeType> m_borders;
    std::vector<SizeType> m_zeros;
    Bvs<SizeType> m_bv;

    KeepLevel() = default;

    KeepLevel(SizeType const size, SizeType const levels, rho_t const& rho) {
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

    SizeType rho(size_t level, size_t i) {
        return (*m_rho)(level, i);
    }

    void set_rho(size_t level, size_t i, SizeType val) {
        // m_rho is already calculated
    }

    std::vector<SizeType>& borders() {
        return m_borders;
    }

    static bool constexpr compute_zeros = is_matrix;
    static bool constexpr compute_rho = false;

    std::vector<SizeType>& zeros() {
        return m_zeros;
    }

    Bvs<SizeType>& bv() {
        return m_bv;
    }
};
