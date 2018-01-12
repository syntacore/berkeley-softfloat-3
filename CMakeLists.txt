cmake_minimum_required(VERSION 3.0)

project(berkeley-softfloat-3)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
# set(CMAKE_C_STANDARD 99)

set(SOFTFLOAT_INLINE "static __inline" CACHE STRING "inline")
set(SOFTFLOAT_INLINE_LEVEL "5" CACHE STRING "INLINE_LEVEL")
set(SOFTFLOAT_FAST_DIV64TO32 OFF CACHE BOOL "SOFTFLOAT_FAST_DIV64TO32")
# set(SOFTFLOAT_FAST_INT64 ON CACHE BOOL "SOFTFLOAT_FAST_INT64")
# set(SOFTFLOAT_LITTLEENDIAN ON CACHE BOOL "LITTLEENDIAN")
# set(SOFTFLOAT_THREAD_LOCAL "" CACHE STRING "thread local")
# set(SOFTFLOAT_FAST_DIV32TO16 ON CACHE BOOL "SOFTFLOAT_FAST_DIV32TO16")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
set(INC_SOFTFLOAT
    ${CMAKE_CURRENT_SOURCE_DIR}/include/softfloat/functions.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/softfloat/types.h
)

add_definitions(-DINLINE=${SOFTFLOAT_INLINE})
# add_definitions(-DTHREAD_LOCAL=${SOFTFLOAT_THREAD_LOCAL})
add_definitions(-DINLINE_LEVEL=${SOFTFLOAT_INLINE_LEVEL})

add_compile_options(-Wall)
if(MSVC) 
    add_compile_options(/wd4820 /wd4514 /MP)
endif()

add_subdirectory(source)
