//
// Created by rgerst on 12.09.18.
//


#pragma once

#include <misaxx-kidney-glomeruli/module_interface.h>
#include <misaxx/core/misa_task.h>

namespace misaxx_kidney_glomeruli {
    struct segmentation3d_base : public misaxx::misa_task {

        misaxx::ome::misa_ome_tiff m_input_segmented2d;
        misaxx::ome::misa_ome_tiff m_output_segmented3d;

        using misaxx::misa_task::misa_task;

    };
}