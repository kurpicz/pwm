################################################################################
# external/wavelet_shun/serialWT/CMakeLists.txt
#
# Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

add_executable(serialWT
  WT.cpp
  WTTime.cpp)

target_link_libraries(serialWT tlx_command_line dl)

if(ENABLE_MALLOC_COUNT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi")
  target_sources(serialWT PUBLIC ${PROJECT_SOURCE_DIR}/src/malloc_count)
  target_compile_definitions(serialWT PUBLIC MALLOC_COUNT)
  target_include_directories(serialWT PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  )
endif(ENABLE_MALLOC_COUNT)

target_include_directories(serialWT PUBLIC
  $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CURRENT_SOURCE_DIR}>
)

################################################################################
