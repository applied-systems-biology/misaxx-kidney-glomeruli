//
// Created by rgerst on 11.09.18.
//


#pragma once

#include <misaxx/core/misa_task.h>
#include <misaxx-kidney-glomeruli/module_interface.h>
#include <misaxx/ome/accessors/misa_ome_plane.h>
#include <misaxx/imaging/coixx/image.h>

namespace misaxx_kidney_glomeruli {
    struct segmentation2d_base : public misaxx::misa_task {

        misaxx::ome::misa_ome_plane<coixx::images::mask> m_input_tissue;
        misaxx::ome::misa_ome_plane<coixx::images::grayscale_float> m_input_autofluoresence;
        misaxx::ome::misa_ome_plane<coixx::images::mask> m_output_segmented2d;

        using misaxx::misa_task::misa_task;

    };
}