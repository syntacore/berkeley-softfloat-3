cmake_minimum_required(VERSION 3.0)

project(berkeley-softfloat-3)

include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(SOFTFLOAT_FAST_DIV64TO32 OFF CACHE BOOL "SOFTFLOAT_FAST_DIV64TO32")

set(SOFTFLOAT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(INC_SOFTFLOAT
    ${CMAKE_CURRENT_SOURCE_DIR}/include/softfloat/functions.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/softfloat/types.h
)

add_library(softfloat-build-defs INTERFACE)

if(IS_BIG_ENDIAN)
    target_compile_definitions(softfloat-build-defs INTERFACE BIG_ENDIAN)
endif()

target_compile_features(softfloat-build-defs INTERFACE
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
target_include_directories(softfloat-build-defs INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(MSVC) 
    target_compile_options(softfloat-build-defs INTERFACE /wd4820 /wd4514 /MP)
else()
    target_compile_options(softfloat-build-defs INTERFACE -Wall)
endif()

add_subdirectory(source)
