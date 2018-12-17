//
// Created by rgerst on 11.09.18.
//


#pragma once

#include <misaxx/workers/misa_task.h>
#include <misaxx_kidney_glomeruli/kidney_glomeruli.h>
#include <misaxx_ome/accessors/misa_ome_plane.h>
#include <coixx/image.h>

namespace misaxx_kidney_glomeruli {
    struct segmentation2d_base : public misaxx::misa_task<kidney_glomeruli> {

        misaxx_ome::misa_ome_plane<coixx::images::mask> m_input_tissue;
        misaxx_ome::misa_ome_plane<coixx::images::grayscale_float> m_input_autofluoresence;
        misaxx_ome::misa_ome_plane<coixx::images::mask> m_output_segmented2d;

        using misaxx::misa_task<kidney_glomeruli>::misa_task;

    };
}