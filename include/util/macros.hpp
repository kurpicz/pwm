/*******************************************************************************
 * include/util/macros.hpp
 *
 * Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#if defined(__GNUC__) || defined(__clang__)
#define PWM_LIKELY(c)   __builtin_expect((c), 1)
#define PWM_UNLIKELY(c) __builtin_expect((c), 0)
#else
#define PWM_LIKELY(c)   c
#define PWM_UNLIKELY(c) c
#endif

#if defined(__GNUC__) || defined(__clang__)
#define PWM_ATTRIBUTE_PACKED __attribute__ ((packed))
#else
#define PWM_ATTRIBUTE_PACKED
#endif

/******************************************************************************/
