cmake_minimum_required(VERSION 3.11)
project(misaxx-kidney-glomeruli VERSION 0.1.0)

find_package(LEMON REQUIRED)
find_package(misaxx-tissue REQUIRED)
find_package(misaxx-ome REQUIRED)

add_library(misaxx-kidney-glomeruli
        include/misaxx-kidney-glomeruli/misaxx_kidney_glomeruli_module.h
        include/misaxx-kidney-glomeruli/kidney_glomeruli.h
        src/misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_base.h
        src/misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_klingberg.h
        src/misaxx-kidney-glomeruli/algorithms/segmentation3d/segmentation3d_klingberg.h
        src/misaxx-kidney-glomeruli/algorithms/segmentation3d/segmentation3d_base.h
        src/misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_local_otsu.h
        src/misaxx-kidney-glomeruli/algorithms/quantification/quantification_base.h
        src/misaxx-kidney-glomeruli/algorithms/quantification/quantification_klingberg.h
        include/misaxx-kidney-glomeruli/attachments/glomerulus.h
        include/misaxx-kidney-glomeruli/attachments/glomeruli.h
        src/misaxx-kidney-glomeruli/algorithms/quantification/quantification_constrained_klingberg.h
        src/misaxx-kidney-glomeruli/kidney_glomeruli.cpp
        src/misaxx-kidney-glomeruli/misaxx_kidney_glomeruli_module.cpp
        src/misaxx-kidney-glomeruli/algorithms/quantification/quantification_constrained_klingberg.cpp
        src/misaxx-kidney-glomeruli/attachments/glomerulus.cpp
        src/misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_klingberg.cpp
        src/misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_local_otsu.cpp
        src/misaxx-kidney-glomeruli/algorithms/quantification/quantification_klingberg.cpp
        src/misaxx-kidney-glomeruli/attachments/glomeruli.cpp
        src/misaxx-kidney-glomeruli/algorithms/segmentation3d/segmentation3d_klingberg.cpp)

target_link_libraries(misaxx-kidney-glomeruli PUBLIC misaxx::misaxx-tissue)
target_link_libraries(misaxx-kidney-glomeruli PUBLIC ${LEMON_LIBRARY})

misaxx_make_shared(misaxx misaxx-kidney-glomeruli)
misaxx_add_default_executable(misaxx-kidney-glomeruli)
misaxx_add_compiler_warnings(misaxx-kidney-glomeruli)
