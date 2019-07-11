cmake_minimum_required(VERSION 3.0)

add_subdirectory(FAST_INT64)
add_subdirectory(NON_FAST_INT64)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(-DSOFTFLOAT_FAST_DIV32TO16)

set(INC_TARGET ${CMAKE_CURRENT_SOURCE_DIR}/Intel8086.hpp)

add_subdirectory(FPU)
add_subdirectory(SSE)
