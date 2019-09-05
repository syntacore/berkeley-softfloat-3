cmake_minimum_required(VERSION 3.8)

add_library(softfloat-${SOFTFLOAT_MODEL}-nfi64 EXCLUDE_FROM_ALL
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128M.cpp
)
target_link_libraries(softfloat-${SOFTFLOAT_MODEL}-nfi64
PUBLIC
	softfloat_interface_nonfast_int64
PRIVATE
	softfloat_x86_nonfast_int64
	softfloat_x86_FPU
)
set_target_properties(softfloat-${SOFTFLOAT_MODEL}-nfi64 PROPERTIES FOLDER "softfloat/${SOFTFLOAT_MODEL}")
