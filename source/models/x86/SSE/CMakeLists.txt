cmake_minimum_required(VERSION 3.0)

add_library(softfloat-x86-SSE EXCLUDE_FROM_ALL)
target_sources(softfloat-x86-SSE
PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
)
target_include_directories(softfloat-x86-SSE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(softfloat-x86-SSE
PUBLIC
	softfloat_interface
PRIVATE
	softfloat_x86
)
set_target_properties(softfloat-x86-SSE PROPERTIES FOLDER "libs/softfloat")
