cmake_minimum_required(VERSION 3.0)

project(berkeley-softfloat-3)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(SOFTFLOAT_FAST_DIV64TO32 OFF CACHE BOOL "SOFTFLOAT_FAST_DIV64TO32")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
set(INC_SOFTFLOAT
    ${CMAKE_CURRENT_SOURCE_DIR}/include/softfloat/functions.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/softfloat/types.h
)

add_compile_options(-Wall)
if(MSVC) 
    add_compile_options(/wd4820 /wd4514 /MP)
endif()

add_subdirectory(source)
