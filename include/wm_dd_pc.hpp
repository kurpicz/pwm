/*******************************************************************************
 * include/wm_dd_pc.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "wx_dd_pc.hpp"

template <typename AlphabetType, typename SizeType = uint64_t>
using wm_dd_pc = wx_dd_pc<AlphabetType, true, SizeType>;
