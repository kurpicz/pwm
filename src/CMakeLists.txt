set(INTERNAL_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/benchmark.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/dd_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/dd_prefix_sorting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/naive.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/parallel_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/parallel_prefix_sorting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/sequential_prefix_counting.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/sequential_prefix_sorting.cpp)

set(INTERNAL_LIBS dl)

add_executable(benchmark ${INTERNAL_SRCS})
target_link_libraries(benchmark ${INTERNAL_LIBS})
