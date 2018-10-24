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

template <typename T>
class span {

public:
  inline span(T* ptr, size_t size) : ptr_(ptr), size_(size) { }

  inline size_t size() const {
    return size_;
  }

  inline T& operator[](size_t i) const {
    DCHECK(i < size());
    return ptr_[i];
  }

  inline T* data() const {
    return ptr_;
  }

  inline span slice(size_t start, size_t end) const {
    DCHECK(start <= m_size);
    DCHECK(end <= m_size);
    DCHECK(start <= end);
    return span{ptr_ + start, end - start};
  }

private:
  T* ptr_;
  size_t size_;
};

/******************************************************************************/
