//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/misa_module_declaration.h>
#include <misaxx/pdata/misa_image_stack.h>
#include <misaxx/pdata/misa_json_file.h>
#include <misaxx/pdata/misa_exportable_meta_data.h>
#include <misaxx_tissue/tissue_module.h>
#include "metadata/glomeruli.h"

using namespace coixx;

namespace misaxx::module::kidney_glomeruli_detection {
    struct kidney_glomeruli : public misa_module_declaration {

        data<misa_image_stack<images::grayscale_float >> m_input_autofluorescence;
        data<misa_image_stack<images::mask>> m_output_segmented2d;
        data<misa_image_stack<images::labels >> m_output_segmented3d;
        data<misa_exportable_meta_data> m_output_quantification;
        misa_submodule<tissue_detection::tissue_detection_module> m_tissue_detection;

        void init_data() override {

            import_from_filesystem(m_input_autofluorescence, "/");
            process(m_output_segmented2d, m_input_autofluorescence, "glomeruli2d");
            process(m_output_segmented3d, m_output_segmented2d, "glomeruli3d");
            export_to_filesystem(m_output_quantification, "quantified.json");

            // Init the submodule
            m_tissue_detection.definition().m_input_autofluorescence = m_input_autofluorescence;
            init_submodule(m_tissue_detection, "tissue");
        }

        glomeruli &get_glomeruli() {
            return m_output_quantification->access<glomeruli>();
        }

        glomerulus &get_glomerulus(int label) {
            return get_glomeruli().data[label];
        }
    };
}