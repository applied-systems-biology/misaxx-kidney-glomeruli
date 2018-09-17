//
// Created by rgerst on 17.09.18.
//


#pragma once

#include "../kidney_glomeruli.h"

namespace misaxx::module::kidney_glomeruli_detection::quantification {
    struct quantification_base : public misa_task<kidney_glomeruli> {

        misa_data<misa_image_stack<images::labels>> m_input_segmented3d;

        using misa_task<kidney_glomeruli>::misa_task;
    };
}