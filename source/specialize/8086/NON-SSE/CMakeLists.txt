cmake_minimum_required(VERSION 3.0)

set(SRC_SPECIALIZE-8086_NON-SSE
    ${SRC_SPECIALIZE-8086}
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF16UI.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF32UI.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF64UI.cpp
)

add_subdirectory(FAST_INT64)
add_subdirectory(NON_FAST_INT64)
