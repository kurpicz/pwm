################################################################################
# external/wavelet_construction/sdslWT/CMakeLists.txt
#
# Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

add_executable(sdslWT
  WT.cpp
  WTTime.cpp)

target_link_libraries(sdslWT tlx_command_line ${SDSL})

# if(ENABLE_MALLOC_COUNT)
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi")
#   target_sources(sdslWT PUBLIC ${PROJECT_SOURCE_DIR}/src/malloc_count)
#   target_compile_definitions(sdslWT PUBLIC MALLOC_COUNT)
#   target_include_directories(sdslWT PUBLIC
#     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
#     $<INSTALL_INTERFACE:${PROJECT_SOURCE_DIR}/include>
#   )
# endif(ENABLE_MALLOC_COUNT)

find_package(SDSL REQUIRED)

target_include_directories(sdslWT PUBLIC
  ${SDSL_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CURRENT_SOURCE_DIR}>
)

################################################################################

