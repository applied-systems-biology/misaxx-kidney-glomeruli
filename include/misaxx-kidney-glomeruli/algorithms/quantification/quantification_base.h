//
// Created by rgerst on 17.09.18.
//


#pragma once

#include <misaxx/misa_task.h>
#include <misaxx-kidney-glomeruli/kidney_glomeruli.h>
#include <misaxx-ome/accessors/misa_ome_tiff.h>

namespace misaxx_kidney_glomeruli {
    struct quantification_base : public misaxx::misa_task {

        misaxx_ome::misa_ome_tiff<coixx::images::labels> m_input_segmented3d;

        using misaxx::misa_task::misa_task;
    };
}