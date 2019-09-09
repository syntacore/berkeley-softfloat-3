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

add_library(softfloat-riscv EXCLUDE_FROM_ALL)
target_link_libraries(softfloat-riscv
PUBLIC
	softfloat_interface
PRIVATE
	softfloat_imp-riscv
	softfloat_public
)
set_target_properties(softfloat-riscv PROPERTIES FOLDER "libs/softfloat")

add_library(softfloat-x86-FPU EXCLUDE_FROM_ALL)
target_link_libraries(softfloat-x86-FPU
PUBLIC
	softfloat_interface
PRIVATE
	softfloat_imp-x86-FPU
	softfloat_public
)
set_target_properties(softfloat-x86-FPU PROPERTIES FOLDER "libs/softfloat")

add_library(softfloat-x86-SSE EXCLUDE_FROM_ALL)
target_link_libraries(softfloat-x86-SSE
PUBLIC
	softfloat_interface
PRIVATE
	softfloat_imp-x86-SSE
	softfloat_public
)
set_target_properties(softfloat-x86-SSE PROPERTIES FOLDER "libs/softfloat")

add_custom_target(softfloat-all
DEPENDS
	softfloat-riscv
	softfloat-x86-FPU
	softfloat-x86-SSE
)
set_target_properties(softfloat-all PROPERTIES FOLDER "libs/softfloat")

if(DOXYGEN_FOUND)
	set(DOXYGEN_STRIP_FROM_PATH ${CMAKE_CURRENT_SOURCE_DIR}/)
	set(DOXYGEN_WARN_LOGFILE ${CMAKE_CURRENT_BINARY_DIR}/doxygen.log)
	set(DOXYGEN_GENERATE_TAGFILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.tag)
	doxygen_add_docs(doxygen-${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR})
endif(DOXYGEN_FOUND)
