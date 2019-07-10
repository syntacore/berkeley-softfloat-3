cmake_minimum_required(VERSION 3.8)

set(SOFTFLOAT_MODEL 8086-SSE)

set(SRC_TARGET_NON_FAST_INT64
    ${SRC_TARGET-8086_NON_FAST_INT64}
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128M.cpp
)

set(ALL_SRC
    ${INC_SOFTFLOAT}
    ${INC_OTHERS}
    ${SRC_OTHERS}
    ${SRC_OTHERS_NON_FAST_INT64}
    ${SRC_PRIMITIVES_NON_FAST_INT64}
    ${SRC_PRIMITIVES}
    ${INC_PRIMITIVES}
    ${INC_TARGET}
    ${SRC_TARGET_NON_FAST_INT64}
)
source_group(TREE ${SOFTFLOAT_ROOT} FILES ${ALL_SRC})

add_library(softfloat-${SOFTFLOAT_MODEL}-nfi64 EXCLUDE_FROM_ALL ${ALL_SRC})
target_link_libraries(softfloat-${SOFTFLOAT_MODEL}-nfi64 PRIVATE softfloat-build-defs)
set_target_properties(softfloat-${SOFTFLOAT_MODEL}-nfi64 PROPERTIES FOLDER "softfloat/${SOFTFLOAT_MODEL}")
