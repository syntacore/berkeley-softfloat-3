cmake_minimum_required(VERSION 3.8)

set(SOFTFLOAT_MODEL RISCV)
add_library(softfloat_riscv INTERFACE)
target_include_directories(softfloat_riscv INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(softfloat_riscv INTERFACE SOFTFLOAT_FAST_DIV32TO16)
target_sources(softfloat_riscv INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/model.hpp)
target_link_libraries(softfloat_riscv INTERFACE softfloat_others)

add_library(softfloat-${SOFTFLOAT_MODEL}-fi64 EXCLUDE_FROM_ALL)
target_link_libraries(softfloat-${SOFTFLOAT_MODEL}-fi64
PUBLIC
	softfloat_interface_fast_int64
PRIVATE
	softfloat_others_fast_int64
	softfloat_riscv
)
set_target_properties(softfloat-${SOFTFLOAT_MODEL}-fi64 PROPERTIES FOLDER "softfloat/${SOFTFLOAT_MODEL}")

add_library(softfloat-${SOFTFLOAT_MODEL}-nfi64 EXCLUDE_FROM_ALL)
target_link_libraries(softfloat-${SOFTFLOAT_MODEL}-nfi64
PUBLIC
	softfloat_interface_nonfast_int64
PRIVATE
	softfloat_others_nonfast_int64
	softfloat_riscv
)
set_target_properties(softfloat-${SOFTFLOAT_MODEL}-nfi64 PROPERTIES FOLDER "softfloat/${SOFTFLOAT_MODEL}")
