//
// Created by rgerst on 17.12.18.
//

#include <misaxx-kidney-glomeruli/algorithms/quantification/quantification_klingberg.h>
#include <misaxx-coixx/objects/label_pixel_count.h>
#include <misaxx-coixx/toolbox/toolbox_objects.h>

using namespace misaxx;
using namespace misaxx_ome;
using namespace misaxx_kidney_glomeruli;
using namespace coixx;

void quantification_klingberg::misa_work() {

    glomeruli result;

    for(const auto &plane : m_input_segmented3d) {
        auto access = plane.access_readonly();

        using namespace coixx::toolbox;

        objects::label_properties<label_pixel_count> component_properties(access.get());

        for(const auto& [group, glom_properties] : component_properties) {

            if(group == 0)
                continue;

            glomerulus glom;
            glom.pixels.count += glom_properties.get<label_pixel_count>().pixels;
            result.data[group] = std::move(glom);
        }
    }

    // Calculate the properties of the glomeruli
    double glomerulus_min_volume = 4.0 / 3.0 * M_PI * pow(m_glomeruli_min_rad, 3);
    double glomerulus_max_volume = 4.0 / 3.0 * M_PI * pow(m_glomeruli_max_rad, 3);

    for(auto &kv : result.data) {
        glomerulus &glom = kv.second;
        glom.label = kv.first;
        glom.volume = glom.pixels.get_volume(module()->m_voxel_size);
        glom.diameter = misa_quantity<double, misa_ome_unit_length<1>>(2 * pow(3.0 / 4.0 * glom.volume.get_value() / M_PI, 1.0 / 3.0),
                                                                       misa_ome_unit_length<1>::ome_unit_type::MICROMETER);
        glom.valid = glom.volume.get_value() >= glomerulus_min_volume && glom.volume.get_value() <= glomerulus_max_volume;
    }

    module()->m_output_quantification.attach(std::move(result));
}