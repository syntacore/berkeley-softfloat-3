cmake_minimum_required(VERSION 3.0)

if(BUILD_SHARED_LIBS)
	add_library(softfloat_imp-x86-SSE STATIC EXCLUDE_FROM_ALL)
	target_sources(softfloat_imp-x86-SSE
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
	)


	target_include_directories(softfloat_imp-x86-SSE PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
	target_link_libraries(softfloat_imp-x86-SSE
	PUBLIC
		softfloat_x86
	PRIVATE
		softfloat_private
	)
	set_target_properties(softfloat_imp-x86-SSE PROPERTIES FOLDER "libs/softfloat")
else()
	add_library(softfloat_imp-x86-SSE INTERFACE)
	target_sources(softfloat_imp-x86-SSE
	INTERFACE
		${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
		${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
	)
	target_include_directories(softfloat_imp-x86-SSE INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
	target_link_libraries(softfloat_imp-x86-SSE
	INTERFACE
		softfloat_x86
		softfloat_private
	)
endif()
