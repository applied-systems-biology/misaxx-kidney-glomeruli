cmake_minimum_required(VERSION 3.11)
project(misaxx-kidney-glomeruli VERSION 0.1.0 DESCRIPTION "MISA++ Kidney Glomeruli Segmentation")

find_package(LEMON REQUIRED)
find_package(misaxx-tissue-segmentation REQUIRED)
find_package(misaxx-imaging-ome REQUIRED)

add_library(misaxx-kidney-glomeruli-segmentation
        include/misaxx-kidney-glomeruli/module.h
        include/misaxx-kidney-glomeruli/module_interface.h
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
        src/misaxx-kidney-glomeruli/module_interface.cpp
        src/misaxx-kidney-glomeruli/module.cpp
        src/misaxx-kidney-glomeruli/algorithms/quantification/quantification_constrained_klingberg.cpp
        src/misaxx-kidney-glomeruli/attachments/glomerulus.cpp
        src/misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_klingberg.cpp
        src/misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_local_otsu.cpp
        src/misaxx-kidney-glomeruli/algorithms/quantification/quantification_klingberg.cpp
        src/misaxx-kidney-glomeruli/attachments/glomeruli.cpp
        src/misaxx-kidney-glomeruli/algorithms/segmentation3d/segmentation3d_klingberg.cpp)

target_link_libraries(misaxx-kidney-glomeruli-segmentation PUBLIC misaxx::misaxx-tissue-segmentation)
target_link_libraries(misaxx-kidney-glomeruli-segmentation PUBLIC ${LEMON_LIBRARY})

set(MISAXX_LIBRARY misaxx-kidney-glomeruli-segmentation)
set(MISAXX_LIBRARY_NAMESPACE misaxx::)
set(MISAXX_API_NAME misaxx_kidney_glomeruli)
set(MISAXX_API_INCLUDE_PATH misaxx-kidney-glomeruli)
set(MISAXX_API_NAMESPACE misaxx_kidney_glomeruli)
misaxx_with_default_module_info()
misaxx_with_default_api()
misaxx_with_default_executable()
