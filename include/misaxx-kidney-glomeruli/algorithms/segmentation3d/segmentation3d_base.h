//
// Created by rgerst on 12.09.18.
//


#pragma once

#include <misaxx-kidney-glomeruli/kidney_glomeruli.h>
#include <misaxx/misa_task.h>
#include <misaxx-coixx/image.h>

namespace misaxx_kidney_glomeruli {
    struct segmentation3d_base : public misaxx::misa_task {

        misaxx_ome::misa_ome_tiff<coixx::images::mask> m_input_segmented2d;
        misaxx_ome::misa_ome_tiff<coixx::images::labels> m_output_segmented3d;

        using misaxx::misa_task::misa_task;

    };
}