cmake_minimum_required(VERSION 3.0)

add_library(softfloat_primitives_fast_int64 INTERFACE)
target_link_libraries(softfloat_primitives_fast_int64
INTERFACE
	softfloat_primitives
	softfloat_interface_fast_int64
)
target_sources(softfloat_primitives_fast_int64
INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/s_mul128To256M.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_mul64To128.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJam128.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJam128Extra.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_sub256M.cpp
)
