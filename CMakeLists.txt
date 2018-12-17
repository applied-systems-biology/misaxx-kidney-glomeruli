cmake_minimum_required(VERSION 3.11)
project(misaxx_kidney_glomeruli)
include("misaxx_cmake/misaxx.cmake")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(LEMON REQUIRED)

add_library(misaxx_kidney_glomeruli
        misaxx_kidney_glomeruli/kidney_glomeruli_detection_module.h
        misaxx_kidney_glomeruli/kidney_glomeruli.h
        misaxx_kidney_glomeruli/algorithms/segmentation2d/segmentation2d_base.h
        misaxx_kidney_glomeruli/algorithms/segmentation2d/segmentation2d_klingberg.h
        misaxx_kidney_glomeruli/algorithms/segmentation3d/segmentation3d_klingberg.h
        misaxx_kidney_glomeruli/algorithms/segmentation3d/segmentation3d_base.h
        misaxx_kidney_glomeruli/algorithms/segmentation2d/segmentation2d_local_otsu.h
        misaxx_kidney_glomeruli/algorithms/quantification/quantification_base.h
        misaxx_kidney_glomeruli/algorithms/quantification/quantification_klingberg.h
        misaxx_kidney_glomeruli/attachments/glomerulus.h
        misaxx_kidney_glomeruli/attachments/glomeruli.h
        misaxx_kidney_glomeruli/algorithms/quantification/quantification_constrained_klingberg.h
        misaxx_kidney_glomeruli/kidney_glomeruli.cpp
        misaxx_kidney_glomeruli/kidney_glomeruli_detection_module.cpp
        misaxx_kidney_glomeruli/algorithms/quantification/quantification_constrained_klingberg.cpp
        misaxx_kidney_glomeruli/attachments/glomerulus.cpp
        misaxx_kidney_glomeruli/algorithms/segmentation2d/segmentation2d_klingberg.cpp
        misaxx_kidney_glomeruli/algorithms/segmentation2d/segmentation2d_local_otsu.cpp
        misaxx_kidney_glomeruli/algorithms/quantification/quantification_klingberg.cpp
        misaxx_kidney_glomeruli/attachments/glomeruli.cpp
        misaxx_kidney_glomeruli/algorithms/segmentation3d/segmentation3d_klingberg.cpp)

target_link_libraries(misaxx_kidney_glomeruli PUBLIC LEMON)

misa_module(misaxx_kidney_glomeruli)
misa_require_submodule(misaxx_kidney_glomeruli misaxx_tissue)
#misa_require_submodule(misaxx_kidney_glomeruli misaxx_ome)