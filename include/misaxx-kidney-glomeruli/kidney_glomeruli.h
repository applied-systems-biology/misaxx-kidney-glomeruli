//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/misa_module_interface.h>
#include <misaxx/accessors/misa_exported_attachments.h>
#include <misaxx-ome/accessors/misa_ome_tiff.h>
#include <misaxx-kidney-glomeruli/attachments/glomeruli.h>
#include <misaxx-tissue/tissue.h>
#include <misaxx-ome/attachments/misa_ome_voxel_size.h>

using namespace coixx;

namespace misaxx_kidney_glomeruli {
    struct kidney_glomeruli : public misaxx::misa_module_interface {

        misaxx_ome::misa_ome_tiff<coixx::images::grayscale_float> m_input_autofluorescence;
        misaxx_ome::misa_ome_tiff<coixx::images::mask> m_output_segmented2d;
        misaxx_ome::misa_ome_tiff<coixx::images::labels > m_output_segmented3d;
        misaxx::misa_exported_attachments m_output_quantification;

        /**
         * Stores the result of the tissue detection
         */
        std::shared_ptr <misaxx_tissue::tissue> m_tissue;

        /**
         * Voxel size obtained from input image
         */
        misaxx_ome::misa_ome_voxel_size m_voxel_size;

        void setup() override;

        glomeruli get_glomeruli();
    };
}