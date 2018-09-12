//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/misa_module_declaration.h>
#include <misaxx/module_data/misa_image_stack.h>
#include <misaxx/module_data/misa_json_file.h>
#include <misaxx/module_data/misa_exportable_meta_data.h>
#include <misaxx_tissue/tissue_module.h>

using namespace coixx;

namespace misaxx::module::kidney_glomeruli_detection {
    struct kidney_glomeruli : public misa_module_declaration {
        data<misa_image_stack<images::grayscale_float >> m_input_autofluorescence = declare_data<misa_image_stack<images::grayscale_float >>("autofluorescence");
        data<misa_image_stack<images::mask>> m_output_segmented2d = declare_data<misa_image_stack<images::mask>>("segmented2d");
        data<misa_image_stack<images::labels >> m_output_segmented3d = declare_data<misa_image_stack<images::labels>>("segmented3d");
        data<misa_exportable_meta_data> m_output_quantification = declare_data<misa_exportable_meta_data>("quantified");
        misa_submodule<tissue_detection::tissue_detection_module> m_tissue_detection = declare_submodule<tissue_detection::tissue_detection_module>("tissue");

        void init_data() override {
            m_input_autofluorescence->from_filesystem(filesystem.imported);
            m_output_quantification->from_filesystem(filesystem.exported / filesystem::as<filesystem::file >("quantified.json"));

            // Init the submodule
            m_tissue_detection.definition().m_input_autofluorescence = m_input_autofluorescence;
            m_tissue_detection.init(*this);

            // Init dependent data
            m_output_segmented2d->from_reference_stack(m_input_autofluorescence);
            m_output_segmented3d->from_reference_stack(m_output_segmented2d);
        }
    };
}