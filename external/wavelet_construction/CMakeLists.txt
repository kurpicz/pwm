################################################################################
# external/wavelet_construction/CMakeLists.txt
#
# Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

add_subdirectory(serialWT)

if (PWM_BUILD_SDSL_WX)
  add_subdirectory(sdslWM)
  add_subdirectory(sdslWT)
endif()

if (PWM_ENABLE_CILK_ALGORITHM)
  add_subdirectory(ddWT)
  #add_subdirectory(ddfsWT)
  add_subdirectory(levelWT)
  add_subdirectory(onlineWT)
  #add_subdirectory(parWT)
  add_subdirectory(recWT)
  add_subdirectory(sortWTOpt)
endif()

################################################################################

