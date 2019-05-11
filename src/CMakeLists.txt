################################################################################
# CMakeLists.txt
#
# Copyright (C) 2016-2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
# Copyright (C) 2017 Marvin Löbel <loebel.marvin@gmail.com>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

set(INTERNAL_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/benchmark.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/dd_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/dd_prefix_counting_single_scan.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/dd_prefix_sorting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/external_parallel_dd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/external_sequential_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/external_parallel_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/external_sequential_prefix_sorting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/naive.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/parallel_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/parallel_prefix_counting_single_scan.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/parallel_prefix_sorting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/sequential_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/sequential_prefix_counting_single_scan.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/sequential_prefix_sorting.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/register/huff_naive.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/huff_dd_naive.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/huff_pc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/huff_dd_pc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/huff_ps.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/huff_dd_ps.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/huff_pc_ss.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/huff_dd_pc_ss.cpp
)

set(INTERNAL_LIBS rt dl)

add_executable(benchmark ${INTERNAL_SRCS})
target_link_libraries(benchmark tlx_command_line stxxl stdc++fs ${INTERNAL_LIBS})

if(ENABLE_MALLOC_COUNT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi")
  target_sources(benchmark PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/malloc_count)
  target_compile_definitions(benchmark PUBLIC MALLOC_COUNT)
endif(ENABLE_MALLOC_COUNT)

target_include_directories(benchmark PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:${PROJECT_SOURCE_DIR}/include/>
)

add_executable(alphabet_count alphabet_count.cpp)
add_executable(fastq_extract fastq_extract.cpp)

################################################################################
