# macro(define_wx prefix cxxdefine)
#     file(GLOB algorithms
#         "${CMAKE_SOURCE_DIR}/include/${prefix}*.hpp"
#     )

#     foreach(loop_var ${algorithms})
#         get_filename_component(WX_NAME "${loop_var}" NAME_WE)
#         #message("${loop_var}")
#         #message("${WX_NAME}")

#         macro(define_wx_variant name definea defineb definec)
#             add_executable(
#                 ${name}_${WX_NAME}
#                 "${CMAKE_SOURCE_DIR}/benchmark/src/construct_${prefix}.cpp"
#             )

#             target_link_libraries(
#                 ${name}_${WX_NAME}
#                 dl
#             )

#             #get_target_property(existing_defs ${WX_NAME} COMPILE_DEFINITIONS)
#             #message("${existing_defs}")

#             set_target_properties(${name}_${WX_NAME} PROPERTIES
#                 COMPILE_DEFINITIONS "${definea};${defineb};${definec}"
#             )
#         endmacro(define_wx_variant)

#         define_wx_variant(construct ${cxxdefine}=${WX_NAME} RUNS=5 TIMING=1)
#         define_wx_variant(check     ${cxxdefine}=${WX_NAME} RUNS=1 CHECK=1)
#         define_wx_variant(memory    ${cxxdefine}=${WX_NAME} RUNS=1 MEMORY=1)

#     endforeach(loop_var)
# endmacro(define_wx)

# define_wx(wm WM_TYPE)
# define_wx(wt WT_TYPE)

# include_directories(.)

set(INTERNAL_SRCS
  benchmark.cpp
  dd_prefix_counting.cpp
  dd_prefix_sorting.cpp
  sequential_prefix_counting.cpp
  sequential_prefix_sorting.cpp)

set(INTERNAL_LIBS dl)

add_executable(benchmark ${INTERNAL_SRCS})
target_link_libraries(benchmark ${INTERNAL_LIBS})