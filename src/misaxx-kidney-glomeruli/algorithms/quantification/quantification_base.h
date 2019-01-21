//
// Created by rgerst on 17.09.18.
//


#pragma once

#include <misaxx/core/misa_task.h>
#include <misaxx-kidney-glomeruli/module_interface.h>
#include <misaxx/ome/accessors/misa_ome_tiff.h>

namespace misaxx_kidney_glomeruli {
    struct quantification_base : public misaxx::misa_task {

        misaxx::ome::misa_ome_tiff<coixx::images::labels> m_input_segmented3d;

        using misaxx::misa_task::misa_task;
    };
}