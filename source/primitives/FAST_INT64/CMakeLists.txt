cmake_minimum_required(VERSION 3.0)

set(SRC_PRIMITIVES_FAST_INT64
    ${CMAKE_CURRENT_SOURCE_DIR}/s_add128.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_add256M.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_mul128To256M.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_mul128By32.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_mul64To128.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJam128.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJam128Extra.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJam256M.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_sub256M.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_mul64ByShifted32To128.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJam64Extra.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftLeft128.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftRight128.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftRightJam128.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftRightJam128Extra.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftRightJam64Extra.cpp
PARENT_SCOPE)
