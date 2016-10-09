cmake_minimum_required(VERSION 3.0)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(INC_OTHERS
    ${CMAKE_CURRENT_SOURCE_DIR}/internals.h
)

add_subdirectory(primitives)
add_subdirectory(others)
add_subdirectory(specialize)
