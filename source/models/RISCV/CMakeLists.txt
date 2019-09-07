cmake_minimum_required(VERSION 3.8)

if(BUILD_SHARED_LIBS)
	add_library(softfloat_imp-riscv STATIC EXCLUDE_FROM_ALL)
	target_sources(softfloat_imp-riscv
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
	)
	target_include_directories(softfloat_imp-riscv PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
	target_link_libraries(softfloat_imp-riscv
	PUBLIC
		softfloat_build
	PRIVATE
		softfloat_private
	)
	set_target_properties(softfloat_imp-riscv PROPERTIES FOLDER "libs/softfloat")
else()
	add_library(softfloat_imp-riscv INTERFACE)
	target_sources(softfloat_imp-riscv
	INTERFACE
		${CMAKE_CURRENT_SOURCE_DIR}/model.hpp
	)
	target_include_directories(softfloat_imp-riscv INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
	target_link_libraries(softfloat_imp-riscv
	INTERFACE
		softfloat_build
		softfloat_private
	)
endif()
