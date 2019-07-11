cmake_minimum_required(VERSION 3.0)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(INC_TARGET ${INC_TARGET} ${CMAKE_CURRENT_SOURCE_DIR}/model.hpp)

set(SRC_TARGET_x86_NON-SSE
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF16UI.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF64UI.cpp
)

set(SOFTFLOAT_MODEL x86-FPU)
add_subdirectory(FAST_INT64)
add_subdirectory(NON_FAST_INT64)
