/*******************************************************************************
 * include/wx_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <vector>

#include "common.hpp"
#include "pc.hpp"

template <typename AlphabetType, bool is_matrix, typename SizeType = uint64_t>
class wx_pc {

public:
    static constexpr bool    is_parallel = false;
    static constexpr bool    is_tree     = !is_matrix;
    static constexpr uint8_t word_width  = sizeof(AlphabetType);

    wx_pc() = default;

    wx_pc(const std::vector<AlphabetType>& text,
         const SizeType size,
         const SizeType levels)
    {
        if(text.size() == 0) { return; }

        auto ctx = LevelSinglePass<SizeType, is_matrix>(size, levels);

        pc(text, size, levels, ctx);

        _zeros = std::move(ctx.zeros());
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
}; // class wx_pc
