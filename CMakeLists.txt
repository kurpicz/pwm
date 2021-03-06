################################################################################
# CMakeLists.txt
#
# Copyright (C) 2016-2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
# Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project (pwm)

# Options
option(PWM_BUILD_COMPETITORS
  "Build other published WT/WM construction algorithms for benchmarks" ON)
option(PWM_ENABLE_CILK_ALGORITHM
  "Enable algorithms using cilk. Onlu if PWM_BUILD_COMPETITORS is ON.\
   Requires gcc version < 8.0.0" ON)
option(PWM_BUILD_SDSL_WX
  "Build the WT construction algorithm depending of the SDSL \
   (requires the SDSL to be installed). Only if PWM_BUILD_COMPETITORS is ON" OFF)

include(ExternalProject)
include(find_or_download_package)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  # using Clang or AppleClang
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "3.4")
    message(FATAL_ERROR "clang version 3.4 or greater required!")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  if(${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS "7.3")
    message(FATAL_ERROR "g++ version 7.3 or greater required!")
  endif()
else()
  message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER_ID} is not supported!")
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} -fopenmp -Wall -Wextra -Wpedantic -fdiagnostics-color=auto")
set(CMAKE_CXX_FLAGS_RELEASE
  "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb -DDEBUG")

# Add include checking
include(header_walker)

# Make malloc count optional
option(ENABLE_MALLOC_COUNT "Enable memory measurement." OFF)
option(USE_GET_RUSAGE_EXTERN "Use getrusage for memory measurement in external/wavelet_construction" ON)

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

if (PWM_BUILD_COMPETITORS)
  # Add competitors
  add_subdirectory(external/wavelet_construction)
endif()

find_or_download_package(STXXL STXXL stxxl)

MESSAGE(STATUS "Built Type: " ${CMAKE_BUILD_TYPE} )

################################################################################
