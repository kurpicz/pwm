ExternalProject_Add(
    glog_external
    GIT_REPOSITORY https://github.com/google/glog.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    UPDATE_COMMAND ""
)
ExternalProject_Get_Property(glog_external source_dir install_dir)

# Only build on demand
set_target_properties(glog_external PROPERTIES EXCLUDE_FROM_ALL TRUE)

file(MAKE_DIRECTORY "${install_dir}/include")

set(${package_found_prefix}_CMAKE_DEP glog_external)
set(${package_found_prefix}_LIBRARIES "${install_dir}/lib/libglog.a")
set(${package_found_prefix}_INCLUDE_DIRS "${install_dir}/include")
