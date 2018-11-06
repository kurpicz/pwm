################################################################################
# CMakeLists.txt
#
# Copyright (C) 2016-2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
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

# Compiler flags
set(CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp -Wall -Wextra -Wpedantic -fdiagnostics-color=auto")
set(CMAKE_CXX_FLAGS_RELEASE
  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb")

# Make malloc count optional
option(ENABLE_MALLOC_COUNT "Enable memory measurement." OFF)

if( NOT EXISTS "${CMAKE_SOURCE_DIR}/external/googletest/.git" )
  execute_process(COMMAND  git submodule update --init
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/external/googletest/"
)
endif()

# Add TLX command line parser
add_library(tlx_command_line
  external/tlx/tlx/cmdline_parser.cpp
  external/tlx/tlx/string/parse_si_iec_units.cpp)

target_include_directories(tlx_command_line PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/tlx/>
  $<INSTALL_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/tlx/>
)

# Add tests
set(BUILD_GMOCK OFF CACHE BOOL "")
set(BUILD_GTEST ON CACHE BOOL "")
add_subdirectory(external/googletest)
add_subdirectory(test)

# Add benchmarks and utility programs
add_subdirectory(src)

# Add competitors
add_subdirectory(external/wavelet_construction)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(ExternalProject)
include(find_or_download_package)
find_or_download_package(STXXL STXXL stxxl)

MESSAGE(STATUS "Built Type: " ${CMAKE_BUILD_TYPE} )

################################################################################
