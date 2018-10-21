################################################################################
# external/wavelet_construction/CMakeLists.txt
#
# Copyright (C) 2018 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

find_package(Cilk REQUIRED)

add_subdirectory(ddWT)
add_subdirectory(levelWT)
add_subdirectory(recWT)
add_subdirectory(sdslWT)
add_subdirectory(serialWT)
add_subdirectory(sortWTOpt)

################################################################################

