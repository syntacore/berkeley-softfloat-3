cmake_minimum_required(VERSION 3.8)

add_library(softfloat-riscv EXCLUDE_FROM_ALL)
target_sources(softfloat-riscv
PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
)
target_include_directories(softfloat-riscv PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(softfloat-riscv PRIVATE SOFTFLOAT_FAST_DIV32TO16)
target_link_libraries(softfloat-riscv
PUBLIC
	softfloat_interface
PRIVATE
	softfloat_others
)
set_target_properties(softfloat-riscv PROPERTIES FOLDER "softfloat")
