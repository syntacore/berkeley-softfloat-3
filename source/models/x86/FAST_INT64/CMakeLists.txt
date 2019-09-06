cmake_minimum_required(VERSION 3.0)

target_sources(softfloat_x86
INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNExtF80UI.cpp
)
