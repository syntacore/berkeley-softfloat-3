cmake_minimum_required(VERSION 3.8)

add_library(softfloat-riscv EXCLUDE_FROM_ALL)
target_sources(softfloat-riscv PRIVATE model.hpp)
target_include_directories(softfloat-riscv PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(softfloat-riscv
PUBLIC
	softfloat_interface
PRIVATE
	softfloat_private
	softfloat_public
)
set_target_properties(softfloat-riscv PROPERTIES FOLDER "libs/softfloat")
