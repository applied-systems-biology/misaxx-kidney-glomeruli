//
// Created by rgerst on 17.09.18.
//


#pragma once

#include "quantification_base.h"
#include "misaxx_kidney_glomeruli/metadata/glomerulus.h"

namespace misaxx::module::kidney_glomeruli_detection::quantification {
    struct klingberg : public quantification_base {

        double  m_glomeruli_min_rad = from_algorithm_json_or<double>("glomeruli-min-rad", 15);
        double  m_glomeruli_max_rad = from_algorithm_json_or<double>("glomeruli-max-rad", 65);
        voxel_size m_voxel_size = from_parameter<voxel_size>();

        using quantification_base::quantification_base;

        void work() {
            for(const auto &layer : *m_input_segmented3d) {

                std::cout << "Layer "<< layer.first << std::endl;
                images::grayscale32s img_components = layer.second->load();

                using namespace coixx::toolbox;

                objects::label_properties<label_pixel_count> component_properties(img_components);

                for(const auto& [group, glom_properties] : component_properties) {

                    if(group == 0)
                        continue;

                    glomerulus &glom = module().get_glomerulus(group);
                    glom.pixels += glom_properties.get<label_pixel_count>().pixels;
                }
            }

            // Calculate the properties of the glomeruli
            double glomerulus_min_volume = 4.0 / 3.0 * M_PI * pow(m_glomeruli_min_rad, 3);
            double glomerulus_max_volume = 4.0 / 3.0 * M_PI * pow(m_glomeruli_max_rad, 3);

            for(auto &kv : module().get_glomeruli().data) {
                glomerulus &glom = kv.second;
                glom.label = kv.first;
                glom.volume = glom.pixels * m_voxel_size.volume();
                glom.diameter = 2 * pow(3.0 / 4.0 * glom.volume / M_PI, 1.0 / 3.0);
                glom.valid = glom.volume >= glomerulus_min_volume && glom.volume <= glomerulus_max_volume;
            }
        }
    };
}