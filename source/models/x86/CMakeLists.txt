cmake_minimum_required(VERSION 3.0)

add_library(softfloat_x86 INTERFACE)
target_sources(softfloat_x86
INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/Intel8086.hpp
)
target_compile_definitions(softfloat_x86 INTERFACE SOFTFLOAT_FAST_DIV32TO16)
target_include_directories(softfloat_x86 INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(softfloat_x86 INTERFACE softfloat_others)

if(SOFTFLOAT_FAST_INT64)
	add_subdirectory(FAST_INT64)
else()
	add_subdirectory(NON_FAST_INT64)
endif(SOFTFLOAT_FAST_INT64)

add_subdirectory(FPU)
add_subdirectory(SSE)
