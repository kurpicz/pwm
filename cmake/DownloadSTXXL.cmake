ExternalProject_Add(
    stxxl_external
    GIT_REPOSITORY https://github.com/stxxl/stxxl
	CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    UPDATE_COMMAND ""
)
ExternalProject_Get_Property(stxxl_external source_dir install_dir)

# Only build on demand
# set_target_properties(stxxl_external PROPERTIES EXCLUDE_FROM_ALL TRUE)

file(MAKE_DIRECTORY "${install_dir}/include")

set(${package_found_prefix}_CMAKE_DEP stxxl_external)
set(${package_found_prefix}_INCLUDE_DIRS "${install_dir}/include")
set(${package_found_prefix}_LIBRARIES
    "${install_dir}/lib/libstxxl.a"
	#"${install_dir}/lib/libstxxl_debug.a"
)

