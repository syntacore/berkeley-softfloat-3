cmake_minimum_required(VERSION 3.0)

add_library(softfloat_build INTERFACE)
target_sources(softfloat_build INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/internals.hpp)
target_include_directories(softfloat_build INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(softfloat_build
INTERFACE
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
if(NOT(MSVC))
	target_compile_options(softfloat_build INTERFACE -Wall)
endif()

add_subdirectory(primitives)
add_subdirectory(others)
add_subdirectory(models/x86)
add_subdirectory(models/RISCV)
