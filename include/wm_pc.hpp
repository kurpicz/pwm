#pragma once

#include "wx_pc.hpp"

template <typename AlphabetType, typename SizeType = uint64_t>
using wm_pc = wx_pc<AlphabetType, true, SizeType>;
