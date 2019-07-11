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

add_library(softfloat-${SOFTFLOAT_MODEL}-fi64 EXCLUDE_FROM_ALL
    ${COMMON}
    ${SRC_OTHERS_FAST_INT64}
    ${SRC_PRIMITIVES_FAST_INT64}
)
target_link_libraries(softfloat-${SOFTFLOAT_MODEL}-fi64 PRIVATE softfloat-build-defs)
set_target_properties(softfloat-${SOFTFLOAT_MODEL}-fi64 PROPERTIES FOLDER "softfloat/${SOFTFLOAT_MODEL}")
target_compile_definitions(softfloat-${SOFTFLOAT_MODEL}-fi64 PRIVATE SOFTFLOAT_FAST_INT64)
if(SOFTFLOAT_FAST_DIV64TO32)
    target_compile_definitions(softfloat-${SOFTFLOAT_MODEL}-fi64 PRIVATE SOFTFLOAT_FAST_DIV64TO32)
endif(SOFTFLOAT_FAST_DIV64TO32)

add_library(softfloat-${SOFTFLOAT_MODEL}-nfi64 EXCLUDE_FROM_ALL
    ${COMMON}
    ${SRC_OTHERS_NON_FAST_INT64}
    ${SRC_PRIMITIVES_NON_FAST_INT64}
)
target_link_libraries(softfloat-${SOFTFLOAT_MODEL}-nfi64 PRIVATE softfloat-build-defs)
set_target_properties(softfloat-${SOFTFLOAT_MODEL}-nfi64 PROPERTIES FOLDER "softfloat/${SOFTFLOAT_MODEL}")

source_group(TREE ${SOFTFLOAT_ROOT} FILES
    ${COMMON}
    ${SRC_OTHERS_FAST_INT64}
    ${SRC_OTHERS_NON_FAST_INT64}
    ${SRC_PRIMITIVES_FAST_INT64}
    ${SRC_PRIMITIVES_NON_FAST_INT64}
)
