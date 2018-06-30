/*******************************************************************************
 * include/util/common.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 * Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <stxxl/vector>
#include "util/type_for_bytes.hpp"

template <typename value_type>
using stxxlvector = typename stxxl::VECTOR_GENERATOR<value_type>::result;

template <int word_width>
using external_vector = stxxlvector<typename type_for_bytes<word_width>::type>;
/******************************************************************************/
