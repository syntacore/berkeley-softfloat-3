cmake_minimum_required(VERSION 3.0)

set(SRC_SPECIALIZE_NON_FAST_INT64
    ${SRC_SPECIALIZE-8086_NON_FAST_INT64}
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128M.c
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
    ${INC_SPECIALIZE}
)
source_group(specialize\\common
FILES
    ${SRC_SPECIALIZE}
)
source_group(specialize\\8086_NON-SSE
FILES
    ${SRC_SPECIALIZE-8086_NON-SSE}
)
source_group(specialize\\NonFastInt64
FILES
    ${SRC_SPECIALIZE_NON_FAST_INT64}
)

add_library(softfloat-8086-nfi64 EXCLUDE_FROM_ALL
    ${INC_SOFTFLOAT}
    ${INC_OTHERS}
    ${SRC_OTHERS}
    ${SRC_OTHERS_NON_FAST_INT64}
    ${SRC_PRIMITIVES_NON_FAST_INT64}
    ${SRC_PRIMITIVES}
    ${INC_PRIMITIVES}
    ${SRC_SPECIALIZE}
    ${SRC_SPECIALIZE-8086_NON-SSE}
    ${INC_SPECIALIZE}
    ${SRC_SPECIALIZE_NON_FAST_INT64}
)

set_target_properties(softfloat-8086-nfi64 PROPERTIES FOLDER "softfloat/8086")
