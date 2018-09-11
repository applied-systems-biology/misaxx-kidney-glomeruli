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
        misa_image_stack<images::grayscale_float > m_input_autofluorescence = data<misa_image_stack<images::grayscale_float >>("autofluorescence");
        misa_image_stack<images::mask> m_output_segmented2d = data<misa_image_stack<images::mask>>("segmented2d");
        misa_image_stack<images::mask> m_output_segmented3d = data<misa_image_stack<images::mask>>("segmented3d");
        misa_exportable_meta_data m_output_quantification = data<misa_exportable_meta_data>("quantified");
        misa_submodule<tissue_detection::tissue_detection_module> m_tissue_detection = imported<tissue_detection::tissue_detection_module>("tissue");
    };
}