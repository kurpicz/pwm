/*******************************************************************************
 * include/arrays/span.hpp
 *
 * Copyright (C) 2018 Marvin Löbel <loebel.marvin@gmail.com>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include <cstring>
#include <iostream>

#include "util/common.hpp"
#include "util/debug_assert.hpp"

template<typename T>
class span {
  T* m_ptr;
  size_t m_size;
public:
  inline span(T* ptr, size_t size):
    m_ptr(ptr), m_size(size) {}
  inline size_t size() const {
    return m_size;
  }
  inline T& operator[](size_t i) const {
    DCHECK(i < size());
    return m_ptr[i];
  }
  inline T* data() const {
    return m_ptr;
  }
  inline span slice(size_t start, size_t end) const {
    DCHECK(start <= m_size);
    DCHECK(end <= m_size);
    DCHECK(start <= end);
    return span { m_ptr + start, end - start };
  }
};


/******************************************************************************/
