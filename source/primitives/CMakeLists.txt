cmake_minimum_required(VERSION 3.0)

add_library(softfloat_primitives INTERFACE)
target_sources(softfloat_primitives
INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/s_addCarryM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_addComplCarryM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_addM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_approxRecip_1Ks.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_approxRecipSqrt_1Ks.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_approxRecipSqrt32_1.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_compare96M.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_compare128M.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_countLeadingZeros8.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_mul64To128.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_mul64To128M.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_mul128MTo256M.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_mul128To256M.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_negXM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_remStepMBy32.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shiftLeftM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJam128.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJam128Extra.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightJamM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shiftRightM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftLeftM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftRightJamM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_shortShiftRightM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_sub1XM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_sub256M.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_subM.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/functions.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/types.hpp
)
target_link_libraries(softfloat_primitives INTERFACE softfloat_build)
