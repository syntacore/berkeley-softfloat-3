cmake_minimum_required(VERSION 3.0)

project(berkeley-softfloat-3)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(SOFTFLOAT_FAST_DIV64TO32 OFF CACHE BOOL "SOFTFLOAT_FAST_DIV64TO32")

set(SOFTFLOAT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(INC_SOFTFLOAT
    ${CMAKE_CURRENT_SOURCE_DIR}/include/softfloat/functions.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/softfloat/types.h
)

add_library(softfloat-build-defs INTERFACE)
target_compile_features(softfloat-build-defs INTERFACE
    cxx_constexpr
)
target_include_directories(softfloat-build-defs INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_compile_options(softfloat-build-defs INTERFACE -Wall)

if(MSVC) 
    target_compile_options(softfloat-build-defs INTERFACE /wd4820 /wd4514 /MP)
endif()

add_subdirectory(source)
