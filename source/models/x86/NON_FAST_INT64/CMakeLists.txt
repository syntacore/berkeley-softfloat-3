cmake_minimum_required(VERSION 3.0)

add_library(softfloat_x86_nonfast_int64 INTERFACE)
target_sources(softfloat_x86_nonfast_int64
INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNExtF80M.cpp
)
target_link_libraries(softfloat_x86_nonfast_int64
INTERFACE
	softfloat_x86
	softfloat_others_nonfast_int64
)
