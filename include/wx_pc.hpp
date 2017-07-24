/*******************************************************************************
 * include/wx_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>

#include "util/common.hpp"
#include "util/pc.hpp"
#include "util/wavelet_structure.hpp"

template <typename AlphabetType, bool is_matrix>
class wx_pc {
    using ctx_t = LevelSinglePass<is_matrix>;

public:
    static constexpr bool    is_parallel = false;
    static constexpr bool    is_tree     = !is_matrix;
    static constexpr uint8_t word_width  = sizeof(AlphabetType);

    wx_pc() = default;

    wx_pc(const std::vector<AlphabetType>& text, const uint64_t size,
        const uint64_t levels) {
        if(text.size() == 0) { return; }

        auto ctx = ctx_t(size, levels);

        pc(text.data(), size, levels, ctx);

        if (ctx_t::compute_zeros)  {
            _zeros = std::move(ctx.zeros());
        }
        _bv = std::move(ctx.bv());
    }

    auto get_bv_and_zeros() const {
        return std::make_pair(_bv.vec(), _zeros);
    }

    auto get_bv() const {
        return _bv.vec();
    }

    wavelet_structure get() && {
        return wavelet_structure(std::move(_bv), std::move(_zeros));
    }
private:
    Bvs _bv;
    std::vector<uint64_t> _zeros;
}; // class wx_pc
