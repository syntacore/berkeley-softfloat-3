cmake_minimum_required(VERSION 3.10)

project(berkeley-softfloat-3)

include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(SOFTFLOAT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(SOFTFLOAT_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(softfloat_interface INTERFACE)
target_sources(softfloat_interface
INTERFACE
	${SOFTFLOAT_INCLUDE}/softfloat/functions.h
	${SOFTFLOAT_INCLUDE}/softfloat/types.h
)
target_include_directories(softfloat_interface INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(softfloat_interface INTERFACE $<$<BOOL:${IS_BIG_ENDIAN}>:BIG_ENDIAN>)

add_library(softfloat_interface_fast_int64 INTERFACE)
target_link_libraries(softfloat_interface_fast_int64 INTERFACE softfloat_interface)
set(SOFTFLOAT_FAST_DIV64TO32 OFF CACHE BOOL "SOFTFLOAT_FAST_DIV64TO32")
target_compile_definitions(softfloat_interface_fast_int64
INTERFACE
	SOFTFLOAT_FAST_INT64
	$<$<BOOL:${SOFTFLOAT_FAST_DIV64TO32}>:SOFTFLOAT_FAST_DIV64TO32>
)

add_library(softfloat_interface_nonfast_int64 INTERFACE)
target_link_libraries(softfloat_interface_nonfast_int64 INTERFACE softfloat_interface)

add_subdirectory(source)
