/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>
#include <type_traits>
#include "common.hpp"

class wavelet_structure {
    Bvs m_bvs;
    std::vector<uint64_t> m_zeros;
public:
    inline wavelet_structure() = default;

    inline wavelet_structure(Bvs&& bvs, std::vector<uint64_t>&& zeros):
        m_bvs(std::move(bvs)), m_zeros(std::move(zeros)) {}

    // Prevent accidental copies
    inline wavelet_structure(wavelet_structure const& other) = delete;
    inline wavelet_structure& operator=(wavelet_structure const& other) = delete;

    // Allow moving
    inline wavelet_structure(wavelet_structure&& other) = default;
    inline wavelet_structure& operator=(wavelet_structure&& other) = default;

    inline std::vector<uint64_t*> const& raw_bvs() const {
        return m_bvs.vec();
    }

    inline std::vector<uint64_t> const& raw_zeros() const {
        return m_zeros;
    }
};

/******************************************************************************/
