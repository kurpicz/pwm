################################################################################
# CMakeLists.txt
#
# Copyright (C) 2016-2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
# Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
project (pwm)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

# Check compiler version
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" AND
   ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "5.4")
  message(FATAL_ERROR "g++ version 5.4 or greater required!")
endif()

set(CXX_STANDARD c++17)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  # using Clang or AppleClang
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "3.4")
    message(FATAL_ERROR "clang version 3.4 or greater required!")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "5.4")
    message(FATAL_ERROR "g++ version 5.4 or greater required!")
  endif()
else()
  message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER_ID} is not supported!")
endif()

# More warnings and debug info
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi")
set(CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp -Wall -Wextra -Wpedantic -fdiagnostics-color=auto")
set(CMAKE_CXX_FLAGS_RELEASE
  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb ")

option(ENABLE_MALLOC_COUNT "Enable memory measurement." OFF)


if( NOT EXISTS "${CMAKE_SOURCE_DIR}/external/googletest/.git" )
  execute_process(COMMAND  git submodule update --init
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/external/googletest/"
)
endif()

# Add include directory to include searches
include_directories(external/tclap/include)
include_directories(include)

set(BUILD_GMOCK OFF CACHE BOOL "")
set(BUILD_GTEST ON CACHE BOOL "")
add_subdirectory(external/googletest)
add_subdirectory(src)
add_subdirectory(test)

MESSAGE(STATUS "Built Type: " ${CMAKE_BUILD_TYPE} )
