#pragma once

#include "wx_dd_pc.hpp"

template <typename AlphabetType, typename SizeType = uint64_t>
using wm_dd_pc = wx_dd_pc<AlphabetType, true, SizeType>;
