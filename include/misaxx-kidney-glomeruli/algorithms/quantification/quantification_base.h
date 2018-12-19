//
// Created by rgerst on 17.09.18.
//


#pragma once

#include <misaxx/workers/misa_task.h>
#include <misaxx-kidney-glomeruli/kidney_glomeruli.h>
#include <misaxx-ome/accessors/misa_ome_tiff.h>
#include <misaxx/workers/misa_task.h>

namespace misaxx_kidney_glomeruli {
    struct quantification_base : public misaxx::misa_task<kidney_glomeruli> {

        misaxx_ome::misa_ome_tiff<coixx::images::labels> m_input_segmented3d;

        using misaxx::misa_task<kidney_glomeruli>::misa_task;
    };
}