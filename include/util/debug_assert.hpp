/*******************************************************************************
 * include/util/debug_assert.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <functional>
#include <iostream>
#include <sstream>

inline void error_with_message(char const* file,
                               uint64_t line,
                               char const* msg,
                               std::function<void(std::ostream&)> extra_msg) {
  std::stringstream ss;
  ss << "In file " << file << ":" << line;
  ss << ":" << std::endl;
  ss << "Check failed: `" << msg << "`" << std::endl;
  std::stringstream ss2;
  extra_msg(ss2);
  auto s = ss2.str();
  if (s.size() > 0) {
    ss << s << std::endl;
  }
  ss << std::endl;
  std::cerr << ss.str() << std::endl;
  std::abort();
}

/// Macro for checking a boolean value
/// and printing a custom error message
#define CHECK_INTERNAL(x, msg, ...)                                            \
  if (!(x)) {                                                                  \
    error_with_message(__FILE__, __LINE__, msg,                                \
                       [&](auto& o) { o << "" __VA_ARGS__; });                 \
  }

// Ensure we have a sane define environment
#if defined(DEBUG) && defined(NDEBUG)
#error Both NDEBUG and DEBUG are set
#endif
#if !defined(DEBUG) && !defined(NDEBUG)
#error Neither NDEBUG nor DEBUG are set
#endif

// Conditionally define DCHECK_INTERNAL
#ifdef DEBUG
#define DCHECK_INTERNAL(x, msg, ...) CHECK_INTERNAL(x, msg, ##__VA_ARGS__)
#else // DEBUG
#define DCHECK_INTERNAL(x, msg, ...)
#endif // DEBUG

// Internal macro for checking a binary operator.
#define DCHECK_OP(op, lhs, rhs)                                                \
  DCHECK_INTERNAL(((lhs) op(rhs)), #lhs " " #op " " #rhs,                      \
                  << "Which is:      " << (lhs) << " " << #op << " " << (rhs))

/// Macro for checking a boolean value. Takes an optional `ostream` message
/// paramter.
///
/// Example:
/// ```
/// DCHECK(predicate(x));
/// DCHECK(y.predicate(), << "invariant in y does not hold");
/// ```
#define DCHECK(x) DCHECK_INTERNAL(x, #x, )
/// Check for equality (==)
#define DCHECK_EQ(x, y) DCHECK_OP(==, x, y)
/// Check for inequality (!=)
#define DCHECK_NE(x, y) DCHECK_OP(!=, x, y)
/// Check for less-than or equal (<=)
#define DCHECK_LE(x, y) DCHECK_OP(<=, x, y)
/// Check for less-than (<)
#define DCHECK_LT(x, y) DCHECK_OP(<, x, y)
/// Check for greater-than or equal (>=)
#define DCHECK_GE(x, y) DCHECK_OP(>=, x, y)
/// Check for greater-than (>)
#define DCHECK_GT(x, y) DCHECK_OP(>, x, y)

/// Macro for checking a boolean value. Takes an optional `ostream` message
/// paramter.
///
/// Example:
/// ```
/// CHECK(predicate(x));
/// CHECK(y.predicate(), << "invariant in y does not hold");
/// ```
#define CHECK(x) CHECK_INTERNAL(x, #x, )

/******************************************************************************/
