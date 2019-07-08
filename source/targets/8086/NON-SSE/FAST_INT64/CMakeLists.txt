cmake_minimum_required(VERSION 3.8)

set(TARGET_NAME 8086)
set(SRC_TARGET_FAST_INT64
    ${SRC_TARGET-8086_FAST_INT64}
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
    ${SRC_TARGET-8086_NON-SSE}
    ${INC_TARGET}
    ${SRC_TARGET_FAST_INT64}
)
source_group(TREE ${SOFTFLOAT_ROOT} FILES ${ALL_SRC})
add_library(softfloat-${TARGET_NAME}-fi64 EXCLUDE_FROM_ALL ${ALL_SRC} )
set_target_properties(softfloat-${TARGET_NAME}-fi64  PROPERTIES FOLDER "softfloat/${TARGET_NAME}"   )
target_compile_definitions(softfloat-${TARGET_NAME}-fi64 PRIVATE SOFTFLOAT_FAST_INT64)
target_include_directories(softfloat-${TARGET_NAME}-fi64 PUBLIC ${SOFTFLOAT_INCLUDE})
if(${SOFTFLOAT_FAST_DIV64TO32})
    target_compile_definitions(softfloat-${TARGET_NAME}-fi64 PRIVATE SOFTFLOAT_FAST_DIV64TO32)
endif(${SOFTFLOAT_FAST_DIV64TO32})
