cmake_minimum_required(VERSION 3.0)

add_library(softfloat_imp-x86-FPU STATIC EXCLUDE_FROM_ALL)
target_sources(softfloat_imp-x86-FPU
PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF16UI.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF64UI.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
)
target_include_directories(softfloat_imp-x86-FPU PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(softfloat_imp-x86-FPU
PUBLIC
	softfloat_x86
PRIVATE
	softfloat_private
)
set_target_properties(softfloat_imp-x86-FPU PROPERTIES FOLDER "libs/softfloat")

add_library(softfloat-x86-FPU EXCLUDE_FROM_ALL)
target_link_libraries(softfloat-x86-FPU
PUBLIC
	softfloat_interface
PRIVATE
	softfloat_imp-x86-FPU
	softfloat_public
)
set_target_properties(softfloat-x86-FPU PROPERTIES FOLDER "libs/softfloat")
