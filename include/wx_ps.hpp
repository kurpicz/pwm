/*******************************************************************************
 * include/wx_ps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>

#include "common.hpp"
#include "ps.hpp"

template <typename AlphabetType, bool is_matrix, typename SizeType = uint64_t>
class wx_ps {
    using ctx_t = LevelSinglePass<SizeType, is_matrix>;

public:
    static constexpr bool    is_parallel = false;
    static constexpr bool    is_tree     = !is_matrix;
    static constexpr uint8_t word_width  = sizeof(AlphabetType);

    wx_ps() = default;

    wx_ps(const std::vector<AlphabetType>& text,
         const SizeType size,
         const SizeType levels)
    {
        if(text.size() == 0) { return; }

        auto ctx = ctx_t(size, levels);

        auto sorted_text = std::vector<AlphabetType>(size);
        ps(text, size, levels, ctx, sorted_text);

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

private:
    Bvs<SizeType> _bv;
    std::vector<SizeType> _zeros;
}; // class wx_ps
