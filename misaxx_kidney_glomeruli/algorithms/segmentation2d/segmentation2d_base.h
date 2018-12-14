//
// Created by rgerst on 11.09.18.
//


#pragma once

#include <misaxx/misa_task.h>
#include "misaxx_kidney_glomeruli/kidney_glomeruli.h"

namespace misaxx_kidney_glomeruli {
    struct segmentation2d_base : public misa_task<kidney_glomeruli> {

        misa_data<misa_image_file<images::mask>> m_input_tissue;
        misa_data<misa_image_file<images::grayscale_float>> m_input_autofluoresence;
        misa_data<misa_image_file<images::mask>> m_output_segmented2d;

        using misa_task::misa_task;

    };
}