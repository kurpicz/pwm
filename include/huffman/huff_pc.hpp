/*******************************************************************************
 * include/huffman/huff_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

template <typename AlphabetType, typename ContextType>
void huff_pc(AlphabetType const* text, const uint64_t size,
  const uint64_t levels, ContextType& ctx) {

  auto& zeros = ctx.zeros();
  auto& borders = ctx.border();
  auto& bv = ctx.bv().vec();
}

/******************************************************************************/
