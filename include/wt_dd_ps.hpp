/*******************************************************************************
 * include/wt_dd_ps.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "wx_dd_ps.hpp"

template <typename AlphabetType, typename SizeType = uint64_t>
using wt_dd_ps = wx_dd_ps<AlphabetType, false, SizeType>;
