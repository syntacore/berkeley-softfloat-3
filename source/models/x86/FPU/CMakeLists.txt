cmake_minimum_required(VERSION 3.0)

add_library(softfloat_x86_FPU EXCLUDE_FROM_ALL)
target_include_directories(softfloat_x86_FPU PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(softfloat_x86_FPU
PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF16UI.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF64UI.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
)
target_link_libraries(softfloat_x86_FPU
PUBLIC
	softfloat_interface
PRIVATE
	softfloat_x86
)
set_target_properties(softfloat_x86_FPU PROPERTIES FOLDER "libs/softfloat")
