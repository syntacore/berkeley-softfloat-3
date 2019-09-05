cmake_minimum_required(VERSION 3.0)

add_library(softfloat_x86_FPU INTERFACE)
target_include_directories(softfloat_x86_FPU INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(softfloat_x86_FPU
INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF16UI.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF64UI.cpp
)
target_link_libraries(softfloat_x86_FPU INTERFACE softfloat_x86)

set(SOFTFLOAT_MODEL x86-FPU)
add_subdirectory(FAST_INT64)
add_subdirectory(NON_FAST_INT64)
