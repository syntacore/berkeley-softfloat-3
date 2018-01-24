cmake_minimum_required(VERSION 3.0)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(INC_TARGET ${INC_TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/target.hpp)
set(SRC_TARGET-8086_SSE ${SRC_TARGET-8086})

add_subdirectory(FAST_INT64)
add_subdirectory(NON_FAST_INT64)
