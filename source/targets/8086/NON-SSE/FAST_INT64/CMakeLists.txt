cmake_minimum_required(VERSION 3.0)

set(SRC_TARGET_FAST_INT64
    ${SRC_TARGET-8086_FAST_INT64}
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
)

source_group(include
FILES
    ${INC_SOFTFLOAT}
)
source_group(primitives
FILES
    ${INC_PRIMITIVES}
)
source_group(primitives\\common
FILES
    ${SRC_PRIMITIVES}
)
source_group(primitives\\FastInt64
FILES
    ${SRC_PRIMITIVES_FAST_INT64}
)
source_group(others
FILES
    ${INC_OTHERS}
)
source_group(others\\common
FILES
    ${SRC_OTHERS}
)
source_group(others\\FastInt64
FILES
    ${SRC_OTHERS_FAST_INT64}
)
source_group(specialize
FILES
    ${INC_TARGET}
)
source_group(specialize\\8086_NON-SSE
FILES
    ${SRC_TARGET-8086_NON-SSE}
)
source_group(specialize\\FastInt64
FILES
    ${SRC_TARGET_FAST_INT64}
)

add_library(softfloat-8086-fi64 EXCLUDE_FROM_ALL
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
set_target_properties(softfloat-8086-fi64  PROPERTIES FOLDER "softfloat/8086"   )

target_compile_definitions(softfloat-8086-fi64 PRIVATE SOFTFLOAT_FAST_INT64)
if(${SOFTFLOAT_FAST_DIV64TO32})
    target_compile_definitions(softfloat-8086-fi64 PRIVATE SOFTFLOAT_FAST_DIV64TO32)
endif(${SOFTFLOAT_FAST_DIV64TO32})
