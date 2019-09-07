cmake_minimum_required(VERSION 3.8)

add_library(softfloat_imp-riscv STATIC EXCLUDE_FROM_ALL)
target_sources(softfloat_imp-riscv
PUBLIC 
	${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
)
target_include_directories(softfloat_imp-riscv PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(softfloat_imp-riscv
PUBLIC
	softfloat_build
PRIVATE
	softfloat_private
)
set_target_properties(softfloat_imp-riscv PROPERTIES FOLDER "libs/softfloat")

add_library(softfloat-riscv EXCLUDE_FROM_ALL)
target_link_libraries(softfloat-riscv
PUBLIC
	softfloat_interface
PRIVATE
	softfloat_imp-riscv
	softfloat_public
)
set_target_properties(softfloat-riscv PROPERTIES FOLDER "libs/softfloat")
