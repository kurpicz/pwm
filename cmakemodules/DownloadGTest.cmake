ExternalProject_Add(
    gtest_external
    GIT_REPOSITORY https://github.com/google/googletest
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    UPDATE_COMMAND ""
)
ExternalProject_Get_Property(gtest_external source_dir install_dir)

# Only build on demand
set_target_properties(gtest_external PROPERTIES EXCLUDE_FROM_ALL TRUE)

# Add thread support
find_package(Threads REQUIRED)
set_target_properties(gtest_external PROPERTIES
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

file(MAKE_DIRECTORY "${install_dir}/include")

set(${package_found_prefix}_CMAKE_DEP gtest_external)
set(${package_found_prefix}_LIBRARIES "${install_dir}/lib/libgtest.a" "${install_dir}/lib/libgtest_main.a")
set(${package_found_prefix}_INCLUDE_DIRS "${install_dir}/include")

