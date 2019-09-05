cmake_minimum_required(VERSION 3.0)

add_library(softfloat_x86_SSE INTERFACE) 
target_include_directories(softfloat_x86_SSE INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(softfloat_x86_SSE INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/model.hpp)
target_link_libraries(softfloat_x86_SSE INTERFACE softfloat_x86)

set(SOFTFLOAT_MODEL SSE)
add_subdirectory(FAST_INT64)
add_subdirectory(NON_FAST_INT64)
