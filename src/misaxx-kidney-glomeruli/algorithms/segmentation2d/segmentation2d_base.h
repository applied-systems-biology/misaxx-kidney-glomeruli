//
// Created by rgerst on 11.09.18.
//


#pragma once

#include <misaxx/misa_task.h>
#include <misaxx-kidney-glomeruli/kidney_glomeruli.h>
#include <misaxx-ome/accessors/misa_ome_plane.h>
#include <misaxx-coixx/image.h>

namespace misaxx_kidney_glomeruli {
    struct segmentation2d_base : public misaxx::misa_task {

        misaxx_ome::misa_ome_plane<coixx::images::mask> m_input_tissue;
        misaxx_ome::misa_ome_plane<coixx::images::grayscale_float> m_input_autofluoresence;
        misaxx_ome::misa_ome_plane<coixx::images::mask> m_output_segmented2d;

        using misaxx::misa_task::misa_task;

    };
}