//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/core/misa_module_interface.h>
#include <misaxx/core/accessors/misa_exported_attachments.h>
#include <misaxx/ome/accessors/misa_ome_tiff.h>
#include <misaxx-kidney-glomeruli/attachments/glomeruli.h>
#include <misaxx-tissue/module_interface.h>
#include <misaxx/ome/attachments/misa_ome_voxel_size.h>
#include <misaxx/imaging/coixx/image.h>

using namespace coixx;

namespace misaxx_kidney_glomeruli {
    struct module_interface : public misaxx::misa_module_interface {

        misaxx::ome::misa_ome_tiff<coixx::images::grayscale_float> m_input_autofluorescence;
        misaxx::ome::misa_ome_tiff<coixx::images::mask> m_output_segmented2d;
        misaxx::ome::misa_ome_tiff<coixx::images::labels > m_output_segmented3d;
        misaxx::misa_exported_attachments m_output_quantification;

        /**
         * Stores the result of the tissue detection
         */
        std::shared_ptr <misaxx_tissue::module_interface> m_tissue;

        /**
         * Voxel size obtained from input image
         */
        misaxx::ome::misa_ome_voxel_size m_voxel_size;

        void setup() override;

        glomeruli get_glomeruli();
    };
}