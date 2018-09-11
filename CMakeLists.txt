cmake_minimum_required(VERSION 3.11)
project(misaxx_kidney_glomeruli)
include(cmake_modules/misa.cmake)
include(cmake_modules/misa_tissue.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(misaxx_kidney_glomeruli misaxx_kidney_glomeruli/kidney_glomeruli_module.h
        misaxx_kidney_glomeruli/kidney_glomeruli.h misaxx_kidney_glomeruli/segmentation2d/segmentation2d_base.h misaxx_kidney_glomeruli/segmentation2d/klingberg.h)
        
target_link_libraries(misaxx_kidney_glomeruli PUBLIC misaxx misaxx_tissue)
target_include_directories(misaxx_kidney_glomeruli INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(misaxx_kidney_glomeruli_exe main.cpp)
target_link_libraries(misaxx_kidney_glomeruli_exe misaxx_kidney_glomeruli)
