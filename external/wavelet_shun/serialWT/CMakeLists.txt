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

target_link_libraries(serialWT tlx_command_line)


target_include_directories(serialWT PUBLIC
  $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CURRENT_SOURCE_DIR}>
)

################################################################################

