include(test/testsuite.cmake)

#list(APPEND BASIC_DEPS dl)
#run_test(vbyte_test     DEPS ${BASIC_DEPS})
#run_test(matrix_tests   DEPS     tudocomp_algorithms BIN_DEPS tudocomp_driver)

macro(define_wx prefix cxxdefine)
    file(GLOB algorithms
        "${CMAKE_SOURCE_DIR}/include/${prefix}*.hpp"
    )

    foreach(loop_var ${algorithms})
        get_filename_component(WX_NAME "${loop_var}" NAME_WE)
        #message("${loop_var}")
        #message("${WX_NAME}")

        if("${WX_NAME}" STREQUAL "wm_cilk_ppc")
            continue()
        endif()
        if("${WX_NAME}" STREQUAL "wt_ppc")
            continue()
        endif()
        if("${WX_NAME}" STREQUAL "wt_pps")
            continue()
        endif()

        set(${cxxdefine} 1)
        configure_file(test_template.cpp ${WX_NAME}_test.cpp)
        run_test(${WX_NAME}_test DEPS dl)
        unset(${cxxdefine})

        #set_target_properties(${WX_NAME}_tests_testrunner PROPERTIES
        #    COMPILE_DEFINITIONS "WX_NAME=${WX_NAME}"
        #)
        #IMPLICIT_DEPENDS_INCLUDE_TRANSFORM "WX_HEADER(%)=<%.hpp>"

    endforeach(loop_var)
endmacro(define_wx)

define_wx(wm WM_TYPE)
define_wx(wt WT_TYPE)

include_directories(.)
