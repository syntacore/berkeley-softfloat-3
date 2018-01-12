cmake_minimum_required(VERSION 3.0)

set(SRC_PRIMITIVES_NON_FAST_INT64
    ${CMAKE_CURRENT_SOURCE_DIR}/s_addCarryM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_addComplCarryM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_addM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_compare128M.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_compare96M.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_mul128MTo256M.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_mul64To128M.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_negXM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_remStepMBy32.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shiftLeftM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJamM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftLeftM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftRightJamM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftRightM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_sub1XM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/s_subM.cpp
PARENT_SCOPE)
