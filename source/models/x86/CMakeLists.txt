cmake_minimum_required(VERSION 3.0)

add_library(softfloat_x86 INTERFACE)
target_sources(softfloat_x86
INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/Intel8086.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNExtF80UI.cpp
)
target_include_directories(softfloat_x86 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(softfloat_x86
INTERFACE
	softfloat_public
	softfloat_private
)

add_subdirectory(FPU)
add_subdirectory(SSE)
