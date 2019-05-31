################################################################################
# external/wavelet_construction/sdslWT/CMakeLists.txt
#
# Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

if (PWM_BUILD_SDSL_WT)

  find_package(SDSL REQUIRED)

  add_executable(sdslWT
    WT.cpp
    WTTime.cpp)

  target_link_libraries(sdslWT tlx_command_line ${SDSL})

  target_include_directories(sdslWT PUBLIC
    ${SDSL_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CURRENT_SOURCE_DIR}>
    )
endif()

################################################################################

