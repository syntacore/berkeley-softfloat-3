cmake_minimum_required(VERSION 3.0)

set(TARGET_NAME 8086-SSE)

set(SRC_TARGET_NON_FAST_INT64
    ${SRC_TARGET-8086_NON_FAST_INT64}
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128M.cpp
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
source_group(primitives\\NonFastInt64
FILES
    ${SRC_PRIMITIVES_NON_FAST_INT64}
)
source_group(others
FILES
    ${INC_OTHERS}
)
source_group(others\\common
FILES
    ${SRC_OTHERS}
)
source_group(others\\NonFastInt64
FILES
    ${SRC_OTHERS_NON_FAST_INT64}
)
source_group(specialize
FILES
    ${INC_TARGET}
)
source_group(specialize\\NonFastInt64
FILES
    ${SRC_TARGET_NON_FAST_INT64}
)

add_library(softfloat-${TARGET_NAME}-nfi64 EXCLUDE_FROM_ALL
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
set_target_properties(softfloat-${TARGET_NAME}-nfi64 PROPERTIES FOLDER "softfloat/${TARGET_NAME}")
