//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/misa_module_declaration.h>
#include <misaxx/accessors/misa_exported_attachments.h>
#include <misaxx_ome/accessors/misa_ome_tiff.h>
#include <misaxx_kidney_glomeruli/attachments/glomeruli.h>
#include <misaxx_tissue/tissue_detection_module.h>
#include <misaxx_ome/attachments/misa_ome_voxel_size.h>

using namespace coixx;

namespace misaxx_kidney_glomeruli {
    struct kidney_glomeruli : public misaxx::misa_module_declaration {

        misaxx_ome::misa_ome_tiff<images::grayscale_float> m_input_autofluorescence;
        misaxx_ome::misa_ome_tiff<images::mask> m_output_segmented2d;
        misaxx_ome::misa_ome_tiff<images::labels > m_output_segmented3d;
        misaxx::misa_exported_attachments m_output_quantification;
        submodule <misaxx_tissue::tissue_detection_module> m_tissue_detection;

        /**
         * Voxel size obtained from input image
         */
        misaxx_ome::misa_ome_voxel_size m_voxel_size;

        void init_data() override;

        glomeruli get_glomeruli();
    };
}