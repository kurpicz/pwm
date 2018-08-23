/*******************************************************************************
 * include/util/debug.hpp
 *
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <iostream>
#include <string>

#include "util/wavelet_structure.hpp"
#include "util/decode.hpp"
#include "huffman/huff_decode.hpp"

[[gnu::unused]] // TODO: C++17 [[maybe_unused]]
static std::string decode_structure(const wavelet_structure& structure) {
    if (structure.is_huffman_shaped()) {
        if (structure.is_tree()) {
            return decode_wt(structure.bvs(), structure.text_size());
        } else {
            return decode_wm(structure.bvs(), structure.zeros(), structure.text_size());
        }
    }  else {
        // TODO
        if (structure.is_tree()) {
            auto& codes = structure.codes<uint8_t, true>();
            return decode_wt_huff<uint8_t>(structure.bvs(), codes);
        } else {
            auto& codes = structure.codes<uint8_t, false>();
            return decode_wm_huff<uint8_t>(structure.bvs(), codes);
        }
        return "";
    }
}

/******************************************************************************/
