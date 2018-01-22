cmake_minimum_required(VERSION 3.0)

add_subdirectory(FAST_INT64)
add_subdirectory(NON_FAST_INT64)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(-DLITTLEENDIAN)
add_definitions(-DSOFTFLOAT_FAST_DIV32TO16)

set(INC_SPECIALIZE ${CMAKE_CURRENT_SOURCE_DIR}/specialize.hpp)
set(SRC_SPECIALIZE-8086
    ${CMAKE_CURRENT_SOURCE_DIR}/s_f16UIToCommonNaN.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_f32UIToCommonNaN.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_f64UIToCommonNaN.cpp
)

add_subdirectory(NON-SSE)
add_subdirectory(SSE)
