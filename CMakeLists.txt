cmake_minimum_required(VERSION 3.10)

project(berkeley-softfloat-3)

include(TestBigEndian)
TEST_BIG_ENDIAN(IS_BIG_ENDIAN)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_library(softfloat_interface INTERFACE)
target_sources(softfloat_interface
INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/include/softfloat/functions.h
	${CMAKE_CURRENT_SOURCE_DIR}/include/softfloat/types.h
)
target_include_directories(softfloat_interface INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(softfloat_interface INTERFACE $<$<BOOL:${IS_BIG_ENDIAN}>:BIG_ENDIAN>)

add_subdirectory(source)
