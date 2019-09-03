################################################################################
# external/wavelet_construction/onlineWT/CMakeLists.txt
#
# Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

set(ONLINE_WT_SRCS
  src/alphabet.cpp
  src/bitarray.cpp
  src/bytearray.cpp
  src/dynarray.cpp
  src/hashset.cpp
  src/strstream.cpp
  src/arrayutil.cpp
  src/bitsandbytes.cpp
  src/csrsbitarray.cpp
  src/dynstr.cpp
  src/huffcode.cpp
  src/wavtree.cpp
  src/binheap.cpp
  src/bitvector.cpp
  src/cstringutil.cpp
  src/hashmap.cpp
  src/mathutil.cpp
)

add_executable(onlineWT
  onlineWT.cpp
  ${ONLINE_WT_SRCS}
)

target_link_libraries(onlineWT tlx_command_line dl "stdc++fs" "c++experimental")

if(ENABLE_MALLOC_COUNT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi")
  target_sources(onlineWT PUBLIC ${PROJECT_SOURCE_DIR}/src/malloc_count)
  target_compile_definitions(onlineWT PUBLIC MALLOC_COUNT)
  target_include_directories(onlineWT PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  )
endif(ENABLE_MALLOC_COUNT)

target_include_directories(onlineWT PUBLIC
  $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CURRENT_SOURCE_DIR}>
)

add_executable(offlineWT
  offlineWT.cpp
  ${ONLINE_WT_SRCS}
)

target_link_libraries(offlineWT tlx_command_line dl "stdc++fs" "c++experimental")

if(ENABLE_MALLOC_COUNT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi")
    target_sources(offlineWT PUBLIC ${PROJECT_SOURCE_DIR}/src/malloc_count)
    target_compile_definitions(offlineWT PUBLIC MALLOC_COUNT)
    target_include_directories(offlineWT PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    )
endif(ENABLE_MALLOC_COUNT)

target_include_directories(offlineWT PUBLIC
  $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CURRENT_SOURCE_DIR}>
)



################################################################################
