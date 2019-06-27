################################################################################
# external/wavelet_construction/sdslWT/CMakeLists.txt
#
# Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

find_package(SDSL REQUIRED)

add_executable(sdslWM
  WT.cpp
  WTTime.cpp)

target_link_libraries(sdslWM tlx_command_line ${SDSL})

target_include_directories(sdslWM PUBLIC
  ${SDSL_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CURRENT_SOURCE_DIR}>
  )

################################################################################

