#include "prefix_counting.hpp"

template <typename AlphabetType, typename SizeType>
using wm_pc = prefix_counting<AlphabetType, SizeType, bit_reverse_permutation>;

template <typename AlphabetType, typename SizeType>
using wt_pc = prefix_counting<AlphabetType, SizeType, identity_function>;