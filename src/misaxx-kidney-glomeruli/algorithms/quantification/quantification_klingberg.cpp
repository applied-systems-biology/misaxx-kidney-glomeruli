//
// Created by rgerst on 17.12.18.
//

#include "quantification_klingberg.h"
#include <misaxx/imaging/coixx/objects/label_pixel_count.h>
#include <misaxx/imaging/coixx/toolbox/toolbox_objects.h>
#include <cmath>

using namespace misaxx;
using namespace misaxx::ome;
using namespace misaxx_kidney_glomeruli;
using namespace coixx;

void quantification_klingberg::work() {
    auto module = get_module_as<module_interface>();

    glomeruli result;
    result.location = misa_location(module->m_output_segmented3d);

    for(const auto &plane : m_input_segmented3d) {
        auto access = plane.access_readonly();

        using namespace coixx::toolbox;

        objects::label_properties<label_pixel_count> component_properties(access.get());

        for(const auto& [group, glom_properties] : component_properties) {

            if(group == 0)
                continue;

            glomerulus glom;
            glom.location = misa_labeled_object_location(module->m_output_segmented3d, group);
            glom.pixels.count += glom_properties.get<label_pixel_count>().pixels;
            result.data[group] = std::move(glom);
        }
    }

    // Calculate the properties of the glomeruli
    double glomerulus_min_volume = 4.0 / 3.0 * M_PI * std::pow(m_glomeruli_min_rad.query(), 3);
    double glomerulus_max_volume = 4.0 / 3.0 * M_PI * std::pow(m_glomeruli_max_rad.query(), 3);

    for(auto &kv : result.data) {
        glomerulus &glom = kv.second;
        glom.label = kv.first;
        glom.volume = glom.pixels.get_volume(module->m_voxel_size);
        glom.diameter = misa_quantity<double, misa_ome_unit_length<1>>(2 * pow(3.0 / 4.0 * glom.volume.get_value() / M_PI, 1.0 / 3.0),
                                                                       misa_ome_unit_length<1>::ome_unit_type::MICROMETER);
        glom.valid = glom.volume.get_value() >= glomerulus_min_volume && glom.volume.get_value() <= glomerulus_max_volume;
    }

    module->m_output_quantification.attach(std::move(result));
}

void quantification_klingberg::create_parameters(misa_parameter_builder &t_parameters) {
    quantification_base::create_parameters(t_parameters);
    m_glomeruli_min_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-min-rad", 15);
    m_glomeruli_max_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-max-rad", 65);
}
