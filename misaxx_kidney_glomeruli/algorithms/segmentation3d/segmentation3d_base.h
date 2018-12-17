//
// Created by rgerst on 12.09.18.
//


#pragma once

#include <misaxx_kidney_glomeruli/kidney_glomeruli.h>
#include <misaxx/workers/misa_task.h>
#include <coixx/image.h>

namespace misaxx_kidney_glomeruli {
    struct segmentation3d_base : public misaxx::misa_task<kidney_glomeruli> {

        misaxx_ome::misa_ome_tiff<coixx::images::mask> m_input_segmented2d;
        misaxx_ome::misa_ome_tiff<coixx::images::labels> m_output_segmented3d;

        using misaxx::misa_task<kidney_glomeruli>::misa_task;

    };
}