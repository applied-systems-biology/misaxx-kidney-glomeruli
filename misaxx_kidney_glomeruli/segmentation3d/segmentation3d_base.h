//
// Created by rgerst on 12.09.18.
//


#pragma once

#include "../kidney_glomeruli.h"

namespace misaxx::module::kidney_glomeruli_detection::segmentation3d {
    struct segmentation3d_base : public misa_task<kidney_glomeruli> {

        misa_data<misa_image_stack<images::mask>> m_input_segmented2d;
        misa_data<misa_image_stack<images::labels>> m_output_segmented3d;

        using misa_task::misa_task;

    };
}