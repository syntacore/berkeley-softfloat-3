cmake_minimum_required(VERSION 3.8)

set(SRC_TARGET_FAST_INT64
    ${SRC_TARGET_x86_FAST_INT64}
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
)

set(ALL_SRC
    ${INC_SOFTFLOAT}
    ${INC_OTHERS}
    ${SRC_OTHERS}
    ${SRC_OTHERS_FAST_INT64}
    ${SRC_PRIMITIVES_FAST_INT64}
    ${SRC_PRIMITIVES}
    ${INC_PRIMITIVES}
    ${SRC_TARGET_x86_NON-SSE}
    ${INC_TARGET}
    ${SRC_TARGET_FAST_INT64}
)
source_group(TREE ${SOFTFLOAT_ROOT} FILES ${ALL_SRC})

add_library(softfloat-${SOFTFLOAT_MODEL}-fi64 EXCLUDE_FROM_ALL ${ALL_SRC} )
target_include_directories(softfloat-${SOFTFLOAT_MODEL}-fi64 PUBLIC ${SOFTFLOAT_INCLUDE})
target_link_libraries(softfloat-${SOFTFLOAT_MODEL}-fi64 PRIVATE softfloat-build-defs)
set_target_properties(softfloat-${SOFTFLOAT_MODEL}-fi64  PROPERTIES FOLDER "softfloat/${SOFTFLOAT_MODEL}"   )
target_compile_definitions(softfloat-${SOFTFLOAT_MODEL}-fi64 PRIVATE SOFTFLOAT_FAST_INT64)
if(SOFTFLOAT_FAST_DIV64TO32)
    target_compile_definitions(softfloat-${SOFTFLOAT_MODEL}-fi64 PRIVATE SOFTFLOAT_FAST_DIV64TO32)
endif(SOFTFLOAT_FAST_DIV64TO32)
