################################################################################
# misc/FindSDSL.cmake
#
# Copyright (C) 2018 Floran Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

# Try to find the SDSL
#
# The following are set after the configuration is done:
# SDSL_FOUND
# SDSL_INCLUDE_DIRS
# SDSL_DIRS
# SDSL

find_path(SDSL_INCLUDE_DIRS sdsl /usr/local/include ~/include ~/local/include)
find_library(SDSL sdsl /usr/local/lib ~/lib ~/local/lib "$ENV{sdsl_ROOT}")

set(SDSL_FOUND TRUE)

if(NOT SDSL)  
  set(SDSL_FOUND FALSE)
  message("-- Finding SDSL failed")
else()  
  get_filename_component(SDSL_DIRS ${SDSL} PATH)
  message("-- Found SDSL Library: ${SDSL_DIRS}")
  message("-- Found SDSL Include: ${SDSL_INCLUDE_DIRS}")
endif(NOT SDSL)

################################################################################

