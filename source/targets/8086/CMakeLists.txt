cmake_minimum_required(VERSION 3.0)

add_subdirectory(FAST_INT64)
add_subdirectory(NON_FAST_INT64)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(-DSOFTFLOAT_FAST_DIV32TO16)

set(INC_TARGET ${CMAKE_CURRENT_SOURCE_DIR}/target.hpp)
set(SRC_TARGET-8086
    ${CMAKE_CURRENT_SOURCE_DIR}/s_f16UIToCommonNaN.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_f32UIToCommonNaN.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_f64UIToCommonNaN.cpp
)

add_subdirectory(NON-SSE)
add_subdirectory(SSE)
