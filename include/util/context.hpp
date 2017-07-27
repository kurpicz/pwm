/*******************************************************************************
 * include/util/context.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "common.hpp"

// TODO: WM/WT abstract that selects zeros and rho
// TODO: flatten vectors where possible, to reduce indirection

// Overwrite information for each level
template<bool is_matrix>
struct LevelSinglePass {
    std::vector<uint64_t> m_hist;
    std::vector<uint64_t> m_borders;
    std::vector<uint64_t> m_zeros;
    Bvs m_bv;
    std::vector<uint64_t> m_bit_reverse;

    LevelSinglePass(uint64_t const size, uint64_t const levels) {
        const auto cur_max_char = (1ull << levels);

        m_hist.reserve(cur_max_char);
        m_hist.resize(cur_max_char, 0);

        if (is_matrix) {
            m_bit_reverse = BitReverse(levels - 1);
        }

        m_borders.reserve(cur_max_char);
        m_borders.resize(cur_max_char, 0);

        m_zeros.reserve(levels);
        m_zeros.resize(levels, 0);

        m_bv = Bvs(size, levels);
    }

    uint64_t hist_size(uint64_t const level) {
        return 1ull << level;
    }

    uint64_t& hist(uint64_t const /*level*/, uint64_t const i) {
        return m_hist[i];
    }

    uint64_t hist(uint64_t const /*level*/, uint64_t const i) const {
        return m_hist[i];
    }

    uint64_t rho(size_t /*level*/, size_t i) {
        if (is_matrix) {
            return m_bit_reverse[i];
        }else {
            return i;
        }
    }

    void set_rho(size_t /*level*/, size_t i, uint64_t val) {
        if (is_matrix) {
            m_bit_reverse[i] = val;
        }
    }

    std::vector<uint64_t>& borders() {
        return m_borders;
    }

    static bool constexpr compute_zeros = is_matrix;
    static bool constexpr compute_rho = is_matrix;

    std::vector<uint64_t>& zeros() {
        return m_zeros;
    }

    Bvs& bv() {
        return m_bv;
    }

    Bvs const& bv() const {
        return m_bv;
    }

    void discard_non_merge_data() {
        // Not used in merge algorithm
    }
};

template <bool is_matrix>
class sliced_single_level_pass {
public:
    sliced_single_level_pass() = default;

    sliced_single_level_pass(uint64_t const size, uint64_t const levels,
        const uint64_t omp_size)
        : hist_(omp_size, std::vector<uint64_t>(1ULL << levels, 0)),
          borders_(omp_size, std::vector<uint64_t>(1ULL << levels, 0)),
          zeros_(levels, 0), bv_(size, levels),
          bit_reverse_(
            is_matrix ? BitReverse(levels - 1) : std::vector<uint64_t>(0)) { }

    uint64_t borders_size(uint64_t const level) {
        return 1ULL << level;
    }

    uint64_t& borders(uint64_t const rank, uint64_t const index) {
        return borders_[rank][index];
    }

    uint64_t borders(uint64_t const rank, uint64_t const index) const {
        return borders_[rank][index];
    }

    uint64_t hist_size(uint64_t const level) {
        return 1ULL << level;
    }

    uint64_t& hist(uint64_t const rank, uint64_t const index) {
        return hist_[rank][index];
    }

    uint64_t hist(uint64_t const rank, uint64_t const index) const {
        return hist_[rank][index];
    }

    uint64_t rho (uint64_t /*level*/, uint64_t const index) {
        if (is_matrix) {
            return bit_reverse_[index];
        } else {
            return index;
        }
    }

    void set_rho(uint64_t /*level*/, uint64_t const index, uint64_t const val) {
        if (is_matrix) {
            bit_reverse_[index] = val;
        }
    }

    static bool constexpr compute_zeros = is_matrix;
    static bool constexpr compute_rho = is_matrix;

    std::vector<uint64_t>& zeros() {
        return zeros_;
    }

    Bvs& bv() {
        return bv_;
    }

    Bvs const& bv() const {
        return bv_;
    }

    void discard_non_merge_data() {
        // Not used in merge algorithms.
    }

private:
    std::vector<std::vector<uint64_t>> hist_;
    std::vector<std::vector<uint64_t>> borders_;

    std::vector<uint64_t> zeros_;
    Bvs bv_;
    std::vector<uint64_t> bit_reverse_;
}; // class sliced_single_level_pass 

/// Keep calculated information for individual levels around
template<bool is_matrix, typename rho_t>
struct KeepLevel {
    std::vector<std::vector<uint64_t>> m_hist;
    rho_t const* m_rho = nullptr;
    std::vector<uint64_t> m_borders;
    std::vector<uint64_t> m_zeros;
    Bvs m_bv;

    KeepLevel() = default;

    KeepLevel(uint64_t const size, uint64_t const levels, rho_t const& rho) {
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

        m_bv = Bvs(size, levels);
    }

    uint64_t hist_size(uint64_t const level) {
        return 1ull << level;
    }

    uint64_t& hist(uint64_t const level, uint64_t const i) {
        return m_hist[level][i];
    }

    uint64_t hist(uint64_t const level, uint64_t const i) const {
        return m_hist[level][i];
    }

    uint64_t rho(size_t level, size_t i) {
        return (*m_rho)(level, i);
    }

    void set_rho(size_t /*level*/, size_t /*i*/, uint64_t /*val*/) {
        // m_rho is already calculated
    }

    std::vector<uint64_t>& borders() {
        return m_borders;
    }

    static bool constexpr compute_zeros = is_matrix;
    static bool constexpr compute_rho = false;

    std::vector<uint64_t>& zeros() {
        return m_zeros;
    }

    Bvs& bv() {
        return m_bv;
    }

    Bvs const& bv() const {
        return m_bv;
    }

    void discard_non_merge_data() {
        drop_me(std::move(m_borders));
    }
};

/******************************************************************************/
