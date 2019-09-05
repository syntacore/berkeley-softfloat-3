cmake_minimum_required(VERSION 3.10)

project(berkeley-softfloat-3)

include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(SOFTFLOAT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(SOFTFLOAT_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(softfloat-interface INTERFACE)
target_include_directories(softfloat-interface INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(IS_BIG_ENDIAN)
    target_compile_definitions(softfloat-interface INTERFACE BIG_ENDIAN)
endif()

add_library(softfloat-interface-fast-int64 INTERFACE)
target_link_libraries(softfloat-interface-fast-int64 INTERFACE softfloat-interface)
target_compile_definitions(softfloat-interface-fast-int64 INTERFACE SOFTFLOAT_FAST_INT64)
set(SOFTFLOAT_FAST_DIV64TO32 OFF CACHE BOOL "SOFTFLOAT_FAST_DIV64TO32")
if(SOFTFLOAT_FAST_DIV64TO32)
    target_compile_definitions(softfloat-interface-fast-int64 INTERFACE SOFTFLOAT_FAST_DIV64TO32)
endif(SOFTFLOAT_FAST_DIV64TO32)

set(INC_SOFTFLOAT
    ${SOFTFLOAT_INCLUDE}/softfloat/functions.h
    ${SOFTFLOAT_INCLUDE}/softfloat/types.h
)

add_library(softfloat-build INTERFACE)
target_compile_features(softfloat-build INTERFACE
    cxx_auto_type
    cxx_alias_templates
    cxx_alignas
    cxx_alignof
    cxx_attributes
    cxx_attribute_deprecated
    cxx_auto_type
    cxx_binary_literals
    cxx_constexpr
    cxx_decltype
    cxx_defaulted_functions
    cxx_delegating_constructors
    cxx_deleted_functions
    cxx_digit_separators
    cxx_enum_forward_declarations
    cxx_explicit_conversions
    cxx_extern_templates
    cxx_final
    cxx_generalized_initializers
    cxx_inline_namespaces
    cxx_lambdas
    cxx_lambda_init_captures
)

if(MSVC) 
    target_compile_options(softfloat-build INTERFACE /wd4820 /wd4514 /MP)
else()
    target_compile_options(softfloat-build INTERFACE -Wall)
endif()

add_subdirectory(source)
