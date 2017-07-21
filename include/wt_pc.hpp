#pragma once

#include "wx_pc.hpp"

template <typename AlphabetType, typename SizeType = uint64_t>
using wt_pc = wx_pc<AlphabetType, false, SizeType>;
