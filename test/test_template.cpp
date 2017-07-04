#include <gtest/gtest.h>

// WX_NAME is defined by the cmake file to be the correct algorithm implementation

//#include "test/util.hpp"

#include "@WX_NAME@.hpp"
#define WX_NAME @WX_NAME@

void foo(WX_NAME<size_t, size_t> bar) {}

TEST(WX_NAME, smoketest) {
}
