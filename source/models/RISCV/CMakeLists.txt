cmake_minimum_required(VERSION 3.8)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(SOFTFLOAT_MODEL RISCV)

add_definitions(-DSOFTFLOAT_FAST_DIV32TO16)

set(INC_TARGET ${CMAKE_CURRENT_SOURCE_DIR}/model.hpp)

set(COMMON
    ${INC_SOFTFLOAT}
    ${INC_OTHERS}
    ${SRC_OTHERS}
    ${INC_PRIMITIVES}
    ${SRC_PRIMITIVES}
    ${INC_TARGET}
)
set(FAST_INT64_ONLY_SRC
    ${SRC_OTHERS_FAST_INT64}
    ${SRC_PRIMITIVES_FAST_INT64}
)
set(NON_FAST_INT64_ONLY_SRC
    ${SRC_OTHERS_NON_FAST_INT64}
    ${SRC_PRIMITIVES_NON_FAST_INT64}
)

source_group(TREE ${SOFTFLOAT_ROOT} FILES
    ${COMMON}
    ${FAST_INT64_ONLY_SRC}
    ${NON_FAST_INT64_ONLY_SRC}
)

add_library(softfloat-${SOFTFLOAT_MODEL}-fi64 EXCLUDE_FROM_ALL
    ${COMMON}
    ${FAST_INT64_ONLY_SRC}
)
target_link_libraries(softfloat-${SOFTFLOAT_MODEL}-fi64
    PUBLIC softfloat-interface-fast-int64
    PRIVATE softfloat-build
)
set_target_properties(softfloat-${SOFTFLOAT_MODEL}-fi64 PROPERTIES FOLDER "softfloat/${SOFTFLOAT_MODEL}")

add_library(softfloat-${SOFTFLOAT_MODEL}-nfi64 EXCLUDE_FROM_ALL
    ${COMMON}
    ${NON_FAST_INT64_ONLY_SRC}
)
target_link_libraries(softfloat-${SOFTFLOAT_MODEL}-nfi64
    PUBLIC softfloat-interface
    PRIVATE softfloat-build
)
set_target_properties(softfloat-${SOFTFLOAT_MODEL}-nfi64 PROPERTIES FOLDER "softfloat/${SOFTFLOAT_MODEL}")
