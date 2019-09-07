cmake_minimum_required(VERSION 3.0)

if(BUILD_SHARED_LIBS)
	add_library(softfloat_imp-x86-FPU STATIC EXCLUDE_FROM_ALL)
	target_sources(softfloat_imp-x86-FPU
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF16UI.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF64UI.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
	)
	target_include_directories(softfloat_imp-x86-FPU PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
	target_link_libraries(softfloat_imp-x86-FPU
	PUBLIC
		softfloat_x86
	PRIVATE
		softfloat_private
	)
	set_target_properties(softfloat_imp-x86-FPU PROPERTIES FOLDER "libs/softfloat")
else()
	add_library(softfloat_imp-x86-FPU INTERFACE)
	target_sources(softfloat_imp-x86-FPU
	INTERFACE
		${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
		${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF16UI.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF64UI.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/s_propagateNaNF128UI.cpp
	)
	target_include_directories(softfloat_imp-x86-FPU INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
	target_link_libraries(softfloat_imp-x86-FPU
	INTERFACE
		softfloat_x86
		softfloat_private
	)
endif()
