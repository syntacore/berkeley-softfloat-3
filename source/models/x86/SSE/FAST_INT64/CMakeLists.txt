cmake_minimum_required(VERSION 3.8)

target_sources(softfloat_x86_SSE
PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
)
