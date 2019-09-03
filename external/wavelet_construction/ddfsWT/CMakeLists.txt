################################################################################
# external/wavelet_construction/ddfsWT/CMakeLists.txt
#
# Copyright (C) 2019 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

set(DDFSWT_SRCS
  bitrank/basic.c
  bitrank/bit_array.c
  bitrank/bitrankw32int.c
  util.c
  dd.c
)

set_source_files_properties(${DDFSWT_SRCS} PROPERTIES LANGUAGE C)

add_executable(ddfsWT
  ${DDFSWT_SRCS})

target_link_libraries(ddfsWT "-std=gnu89 -lrt -fcilkplus -lcilkrts -lm" dl)

if(ENABLE_MALLOC_COUNT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi")
    target_sources(ddfsWT PUBLIC ${PROJECT_SOURCE_DIR}/src/malloc_count)
    target_compile_definitions(ddfsWT PUBLIC MALLOC_COUNT)
    target_include_directories(ddfsWT PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    )
endif(ENABLE_MALLOC_COUNT)

target_include_directories(ddfsWT PUBLIC
  $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CURRENT_SOURCE_DIR}>
)

target_compile_options(ddfsWT PRIVATE -fcilkplus)

################################################################################
