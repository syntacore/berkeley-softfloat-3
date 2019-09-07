cmake_minimum_required(VERSION 3.0)

add_library(softfloat_primitives INTERFACE)
target_sources(softfloat_primitives
INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/s_approxRecip_1Ks.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_approxRecipSqrt_1Ks.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_approxRecipSqrt32_1.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_countLeadingZeros8.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_mul64To128.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_mul128To256M.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJam128.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJam128Extra.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_sub256M.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/functions.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/types.hpp
)
target_link_libraries(softfloat_primitives INTERFACE softfloat_build)

if(NOT(SOFTFLOAT_FAST_INT64))
	add_subdirectory(NON_FAST_INT64)
endif()
