# Custom test target to run the googletest tests
add_custom_target(check)
add_custom_command(
    TARGET check
    POST_BUILD
    COMMENT "All tests were successful!" VERBATIM
)

# Custom test target to just build the googletest tests
add_custom_target(build_check)
add_custom_command(
    TARGET build_check
    POST_BUILD
    COMMENT "All test builds were successful!" VERBATIM
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/stamps)

# will compile and run ${test_target}.cpp
# and add all further arguments as dependencies
macro(generic_run_test test_target test_file
      driver driver_dep register_target register_build_target kind_name construction_test)
    set(options)
    set(oneValueArgs)
    set(multiValueArgs DEPS BIN_DEPS)
    cmake_parse_arguments(TEST_TARGET "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if (${construction_test})
        set(CONSTRUCTION_ALGORITHMS
            ${PROJECT_SOURCE_DIR}/src/dd_prefix_counting.cpp
            ${PROJECT_SOURCE_DIR}/src/dd_prefix_sorting.cpp
            ${PROJECT_SOURCE_DIR}/src/naive.cpp
            ${PROJECT_SOURCE_DIR}/src/parallel_prefix_counting.cpp
            # ${PROJECT_SOURCE_DIR}/src/parallel_prefix_sorting.cpp
            ${PROJECT_SOURCE_DIR}/src/sequential_prefix_counting.cpp
            ${PROJECT_SOURCE_DIR}/src/sequential_prefix_sorting.cpp)
    endif()

    add_executable(${test_target}_testrunner
        EXCLUDE_FROM_ALL
        ${driver}
        ${test_file}
        ${CONSTRUCTION_ALGORITHMS}
    )

    target_link_libraries(${test_target}_testrunner
        ${driver_dep}
        ${TEST_TARGET_DEPS}
    )

    # Runs the test and generates a stamp file on success.
    add_custom_command(
        OUTPUT
            stamps/${test_target}_testrunner.stamp
        DEPENDS
            ${test_target}_testrunner
        COMMAND
            ${test_target}_testrunner
        COMMAND
            cmake -E touch ${CMAKE_CURRENT_BINARY_DIR}/stamps/${test_target}_testrunner.stamp
        WORKING_DIRECTORY
            "${CMAKE_BINARY_DIR}"
        COMMENT
            "Running ${kind_name} ${test_target} ..."
        VERBATIM
    )

    # The test target. Depends on the stamp file to ensure the
    # test is only run if the source changed
    add_custom_target(
        ${test_target}
        DEPENDS
            stamps/${test_target}_testrunner.stamp
    )

    # Hook into check target
    add_custom_command(
        TARGET ${register_target}
        PRE_BUILD
        COMMAND cmake --build . --target ${test_target}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        COMMENT "${kind_name} ${test_target}" VERBATIM
    )

    # Hook into build_check target
    add_custom_command(
        TARGET ${register_build_target}
        PRE_BUILD
        COMMAND cmake --build . --target ${test_target}_testrunner
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        COMMENT "Building ${kind_name} ${test_target}" VERBATIM
    )

    # Ensure binary deps of the testrunner are compiled first
    foreach(bin_dep ${TEST_TARGET_BIN_DEPS})
        add_custom_command(
            TARGET ${test_target}_testrunner
            PRE_BUILD
            COMMAND cmake --build . --target ${bin_dep}
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        )
    endforeach(bin_dep)
endmacro()

macro(run_test test_target construction_test)
generic_run_test(
    ${test_target}
    "${test_target}.cpp"
    "test/test_driver.cpp"
    gtest
    check
    build_check
    "Test"
    "${construction_test}"
    ${ARGN}
)
endmacro()
