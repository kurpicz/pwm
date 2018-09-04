set(INTERNAL_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/benchmark.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/dd_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/dd_prefix_counting_single_scan.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/dd_prefix_sorting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/naive.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/parallel_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/parallel_prefix_counting_single_scan.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/parallel_prefix_sorting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/sequential_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/sequential_prefix_counting_single_scan.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/register/sequential_prefix_sorting.cpp)

if(ENABLE_MALLOC_COUNT)
  add_definitions(-DMALLOC_COUNT)
  list(APPEND INTERNAL_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/malloc_count)
endif(ENABLE_MALLOC_COUNT)

set(INTERNAL_LIBS rt dl)

add_executable(benchmark ${INTERNAL_SRCS})
target_link_libraries(benchmark ${INTERNAL_LIBS})

add_executable(alphabet_count alphabet_count.cpp)
add_executable(fastq_extract fastq_extract.cpp)
