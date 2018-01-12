cmake_minimum_required(VERSION 3.0)

set(SRC_SPECIALIZE
    ${CMAKE_CURRENT_SOURCE_DIR}/extF80M_isSignalingNaN.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/f128M_isSignalingNaN.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/softfloat_raiseFlags.cpp
)

add_subdirectory(8086)
add_subdirectory(RISCV)
