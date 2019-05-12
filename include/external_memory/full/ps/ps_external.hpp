/*******************************************************************************
 * include/util/ps_external.hpp
 *
 * Copyright (C) 2018 Jonas Ellert <jonas.ellert@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "external_memory/full/ps/packed_io.hpp"
#include "external_memory/wavelet_structure_external.hpp"
#include "util/debug_assert.hpp"

#define PSE_VERBOSE if constexpr (false)
#define PSE_BUILD_FUNC static void build(const InputType& text, \
    wavelet_structure_external& result, \
    stats_type& stats);

template <typename InputType, typename stats_type,
          bool is_tree, int word_packing_mode = 2>
struct wx_ps_fe_builder;

template <typename InputType, typename stats_type, int word_packing_mode>
struct wx_ps_fe_builder<InputType, stats_type, false, word_packing_mode> {
  PSE_BUILD_FUNC };

template <typename InputType, typename stats_type, int word_packing_mode>
struct wx_ps_fe_builder<InputType, stats_type, true, word_packing_mode> {
  PSE_BUILD_FUNC };

template <typename InputType, typename stats_type>
struct wx_ps_fe_builder<InputType, stats_type, false, 0> {
  PSE_BUILD_FUNC };

template <typename InputType, typename stats_type>
struct wx_ps_fe_builder<InputType, stats_type, true, 0> {
  PSE_BUILD_FUNC };


#include "external_memory/full/ps/ps_external_packed.hpp"
#include "external_memory/full/ps/ps_external_unpacked.hpp"

/******************************************************************************/
