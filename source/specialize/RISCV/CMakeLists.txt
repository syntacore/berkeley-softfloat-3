cmake_minimum_required(VERSION 3.0)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(-DLITTLEENDIAN)
add_definitions(-DSOFTFLOAT_FAST_DIV32TO16)

set(INC_SPECIALIZE ${CMAKE_CURRENT_SOURCE_DIR}/specialize.h)
set(SRC_SPECIALIZE_RISCV
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF16UI.c
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF32UI.c
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF64UI.c
)
set(SRC_SPECIALIZE_FAST_INT64
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNExtF80UI.c
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.c
)
set(SRC_SPECIALIZE_NON_FAST_INT64
    ${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNExtF80M.c
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
    ${INC_SPECIALIZE}
)
source_group(specialize\\common
FILES
    ${SRC_SPECIALIZE}
)
source_group(specialize\\RISCV
FILES
    ${SRC_SPECIALIZE_RISCV}
)
source_group(specialize\\FastInt64
FILES
    ${SRC_SPECIALIZE_FAST_INT64}
)
source_group(specialize\\NonFastInt64
FILES
    ${SRC_SPECIALIZE_NON_FAST_INT64}
)

add_library(softfloat-RISCV-fi64 EXCLUDE_FROM_ALL
    ${INC_SOFTFLOAT}
    ${INC_OTHERS}
    ${SRC_OTHERS}
    ${SRC_OTHERS_FAST_INT64}
    ${SRC_PRIMITIVES_FAST_INT64}
    ${SRC_PRIMITIVES}
    ${INC_PRIMITIVES}
    ${SRC_SPECIALIZE}
    ${INC_SPECIALIZE}
    ${SRC_SPECIALIZE_FAST_INT64}
    ${SRC_SPECIALIZE_RISCV}
)
add_library(softfloat-RISCV-nfi64 EXCLUDE_FROM_ALL
    ${INC_SOFTFLOAT}
    ${INC_OTHERS}
    ${SRC_OTHERS}
    ${SRC_OTHERS_NON_FAST_INT64}
    ${SRC_PRIMITIVES_NON_FAST_INT64}
    ${SRC_PRIMITIVES}
    ${INC_PRIMITIVES}
    ${SRC_SPECIALIZE}
    ${INC_SPECIALIZE}
    ${SRC_SPECIALIZE_NON_FAST_INT64}
    ${SRC_SPECIALIZE_RISCV}
)

set_target_properties(softfloat-RISCV-fi64 PROPERTIES FOLDER "softfloat/RISCV")
set_target_properties(softfloat-RISCV-nfi64 PROPERTIES FOLDER "softfloat/RISCV")

target_compile_definitions(softfloat-RISCV-fi64 PRIVATE SOFTFLOAT_FAST_INT64)
if(${SOFTFLOAT_FAST_DIV64TO32})
    target_compile_definitions(softfloat-RISCV-fi64 PRIVATE SOFTFLOAT_FAST_DIV64TO32)
endif(${SOFTFLOAT_FAST_DIV64TO32})
