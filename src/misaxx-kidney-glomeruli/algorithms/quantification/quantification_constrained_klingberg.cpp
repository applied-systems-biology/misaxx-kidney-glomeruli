//
// Created by rgerst on 14.12.18.
//

#include <misaxx-kidney-glomeruli/algorithms/quantification/quantification_constrained_klingberg.h>
#include <misaxx-coixx/objects/label_pixel_count.h>
#include <misaxx-coixx/toolbox/toolbox_objects.h>

using namespace misaxx;
using namespace misaxx_ome;
using namespace misaxx_kidney_glomeruli;
using namespace coixx;

void quantification_constrained_klingberg::misa_work() {

    glomeruli result;

    for(size_t layer_index = 0; layer_index < m_input_segmented3d.size(); ++layer_index) {

        images::grayscale32s img_components = m_input_segmented3d.at(layer_index).clone();

        using namespace coixx::toolbox;
        objects::label_properties<label_pixel_count, label_min_max_position> component_properties(img_components);

        for(const auto& [group, glom_properties] : component_properties) {

            if(group == 0)
                continue;

            glomerulus &glom = result.data[group];
            glom.pixels.count += glom_properties.get<label_pixel_count>().pixels;

            auto &bb = glom.bounds;
            const auto &vs = module()->m_voxel_size;
            const auto &lbl_minmax = glom_properties.get<label_min_max_position>();

            // Count z-layer bounding
            bb.include_z(vs.get_size_z().make(layer_index));

            // Count the bounding
            bb.include_x(vs.get_size_x().make(lbl_minmax.min_x));
            bb.include_x(vs.get_size_x().make(lbl_minmax.max_x));
            bb.include_y(vs.get_size_y().make(lbl_minmax.min_y));
            bb.include_y(vs.get_size_y().make(lbl_minmax.max_y));
        }
    }

    // Calculate the properties of the glomeruli

    const double glomerulus_min_volume = 4.0 / 3.0 * M_PI * pow(m_glomeruli_min_rad, 3);
    const double glomerulus_max_volume = 4.0 / 3.0 * M_PI * pow(m_glomeruli_max_rad, 3);
    const double glomerulus_min_diameter = 2 * m_glomeruli_min_rad;
    const double glomerulus_max_diameter = 2 * m_glomeruli_max_rad;

    for(auto &kv : result.data) {

        glomerulus &glom = kv.second;

        glom.label = kv.first;
        glom.volume = glom.pixels.get_volume(module()->m_voxel_size);
        glom.diameter = misa_quantity<double, misa_ome_unit_length<1>>(2 * std::pow(3.0 / 4.0 * glom.volume.get_value() / M_PI, 1.0 / 3.0),
                                                                       misa_ome_unit_length<1>::ome_unit_type::MICROMETER);

        const auto &bb = glom.bounds;

        if(bb.is_valid()) {
            const double feret_x = bb.get_size_x().get_value();
            const double feret_y = bb.get_size_y().get_value();
            const double feret_z = bb.get_size_z().get_value();

            glom.valid = glom.volume.get_value() >= glomerulus_min_volume && glom.volume.get_value() <= glomerulus_max_volume &&
                         bb.is_valid() &&
                         feret_x >= glomerulus_min_diameter && feret_x <= glomerulus_max_diameter &&
                         feret_y >= glomerulus_min_diameter && feret_y <= glomerulus_max_diameter &&
                         feret_z >= glomerulus_min_diameter && feret_z <= glomerulus_max_diameter;
        }
        else {
            glom.valid = false;
        }
    }

    module()->m_output_quantification.attach(std::move(result));
}
