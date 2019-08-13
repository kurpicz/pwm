
################################################################################
# external/wavelet_construction/parWT/CMakeLists.txt
#
# Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

add_executable(parWT
  bitrank/basic.c
  bitrank/bit_array.c
  bitrank/bitrankw32int.c
  pwt.c
  util.c)

target_link_libraries(parWT tlx_command_line cilkrts dl)

if(ENABLE_MALLOC_COUNT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi")
    target_sources(parWT PUBLIC ${PROJECT_SOURCE_DIR}/src/malloc_count)
    target_compile_definitions(parWT PUBLIC MALLOC_COUNT)
    target_include_directories(parWT PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    )
endif(ENABLE_MALLOC_COUNT)

target_include_directories(parWT PUBLIC
  $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CURRENT_SOURCE_DIR}>
)

target_compile_options(parWT PRIVATE -fcilkplus)

################################################################################
