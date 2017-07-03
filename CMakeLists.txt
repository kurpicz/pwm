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

set(CXX_STANDARD c++14)

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

include(ExternalProject)

# More warnings and debug info
# TODO: Add these, they break building glog -Wextra -Wpedantic
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp -fcilkplus -ffast-math -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb ")

# Find Python3
#set(Python_ADDITIONAL_VERSIONS 3)
#find_package(PythonInterp REQUIRED)

# Global dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmakemodules")

# Finders
include(find_or_download_package)

# Usage: find_or_download_package(PACKAGE_NAME, FOUND_FLAG_PREFIX, TARGET_PREFIX)
find_or_download_package(Glog GLOG glog)
find_or_download_package(GTest GTEST gtest)

# Add include directory to include searches
include_directories(include)

add_subdirectory(src)
add_subdirectory(test)

MESSAGE( STATUS "Built Type: " ${CMAKE_BUILD_TYPE} )
