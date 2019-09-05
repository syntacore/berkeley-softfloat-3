cmake_minimum_required(VERSION 3.8)

add_library(softfloat-${SOFTFLOAT_MODEL}-fi64 EXCLUDE_FROM_ALL
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
)
target_link_libraries(softfloat-${SOFTFLOAT_MODEL}-fi64
PUBLIC
	softfloat_interface_fast_int64
PRIVATE
	softfloat_x86_fast_int64
	softfloat_x86_SSE
)
set_target_properties(softfloat-${SOFTFLOAT_MODEL}-fi64 PROPERTIES FOLDER "softfloat/${SOFTFLOAT_MODEL}")
