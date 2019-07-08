cmake_minimum_required(VERSION 3.0)

project(berkeley-softfloat-3)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(SOFTFLOAT_FAST_DIV64TO32 OFF CACHE BOOL "SOFTFLOAT_FAST_DIV64TO32")

set(SOFTFLOAT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(SOFTFLOAT_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(INC_SOFTFLOAT
    ${SOFTFLOAT_INCLUDE}/softfloat/functions.h
    ${SOFTFLOAT_INCLUDE}/softfloat/types.h
)

add_compile_options(-Wall)
if(MSVC) 
    add_compile_options(/wd4820 /wd4514 /MP)
endif()

add_subdirectory(source)
