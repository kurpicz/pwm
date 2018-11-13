find_program(iwyu_path NAMES include-what-you-use iwyu)
if(NOT iwyu_path)
  message(STATUS "Could not find the program include-what-you-use")
else()
  message(STATUS "Found the program include-what-you-use at ${iwyu_path}")
endif()

set(IWYU_CMAKE_DIR ${CMAKE_CURRENT_LIST_DIR})
macro(add_iwyu_diagnostic_for_ccx_target target)
if(iwyu_path)
    set_property(TARGET ${target} PROPERTY CXX_INCLUDE_WHAT_YOU_USE "${iwyu_path};-Xiwyu;--mapping_file=${IWYU_CMAKE_DIR}/iwyu.imp")
endif()
endmacro()
