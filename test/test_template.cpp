#include <gtest/gtest.h>
#include "test/util.hpp"

// WX_NAME is defined by the cmake file to be one of the algorithms
#include "@WX_NAME@.hpp"
#define WX_NAME @WX_NAME@

void foo(WX_NAME<size_t, size_t> bar) {}

TEST(WX_NAME, smoketest) {
    //std::cout << "foo\n";
}
