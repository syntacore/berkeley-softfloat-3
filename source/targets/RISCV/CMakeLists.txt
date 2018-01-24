cmake_minimum_required(VERSION 3.0)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(TARGET_NAME RISCV)

add_definitions(-DSOFTFLOAT_FAST_DIV32TO16)

set(INC_TARGET ${CMAKE_CURRENT_SOURCE_DIR}/target.hpp)
set(SRC_TARGET)
set(SRC_TARGET_NON_FAST_INT64)
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
source_group(others\\FastInt64
FILES
    ${SRC_OTHERS_FAST_INT64}
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

add_library(softfloat-${TARGET_NAME}-fi64 EXCLUDE_FROM_ALL
    ${INC_SOFTFLOAT}
    ${INC_OTHERS}
    ${SRC_OTHERS}
    ${SRC_OTHERS_FAST_INT64}
    ${SRC_PRIMITIVES_FAST_INT64}
    ${SRC_PRIMITIVES}
    ${INC_PRIMITIVES}
    ${INC_TARGET}
    ${SRC_TARGET}
)
set_target_properties(softfloat-${TARGET_NAME}-fi64 PROPERTIES FOLDER "softfloat/${TARGET_NAME}")

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
    ${SRC_TARGET}
)
set_target_properties(softfloat-${TARGET_NAME}-nfi64 PROPERTIES FOLDER "softfloat/${TARGET_NAME}")

target_compile_definitions(softfloat-${TARGET_NAME}-fi64 PRIVATE SOFTFLOAT_FAST_INT64)
if(${SOFTFLOAT_FAST_DIV64TO32})
    target_compile_definitions(softfloat-${TARGET_NAME}-fi64 PRIVATE SOFTFLOAT_FAST_DIV64TO32)
endif(${SOFTFLOAT_FAST_DIV64TO32})
