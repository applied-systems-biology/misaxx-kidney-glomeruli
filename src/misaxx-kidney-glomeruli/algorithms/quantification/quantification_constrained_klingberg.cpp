//
// Created by rgerst on 14.12.18.
//

#include "quantification_constrained_klingberg.h"
#include <cmath>
#include <cv-toolbox/ReadableBMatTypes.h>
#include <cv-toolbox/label_properties.h>

using namespace misaxx;
using namespace misaxx::ome;
using namespace misaxx_kidney_glomeruli;

namespace {
    
    /**
     * Properties collected for labeling
     */
    struct cc_properties {
        size_t pixels = 0;
        int min_x = std::numeric_limits<int>::max();
        int max_x = std::numeric_limits<int>::lowest();
        int min_y = std::numeric_limits<int>::max();
        int max_y = std::numeric_limits<int>::lowest();
        
        void update(int x, int y, int label) {
            ++pixels;
            min_x = std::min(min_x, x);
            min_y = std::min(min_y, y);
            max_x = std::max(max_x, x);
            max_y = std::max(max_y, y);
        }
    };
    
}

void quantification_constrained_klingberg::work() {
    auto module = get_module_as<module_interface>();

    glomeruli result;

    for(size_t layer_index = 0; layer_index < m_input_segmented3d.size(); ++layer_index) {

        cv::images::labels img_components { m_input_segmented3d.at(layer_index).clone() };
        cv::label_properties<cc_properties> component_properties(img_components);

        for(const auto& [group, glom_properties] : component_properties.rows) {

            if(group == 0)
                continue;

            glomerulus &glom = result.data[group];
            glom.pixels.count += glom_properties.pixels;

            auto &bb = glom.bounds;
            const auto &vs = module->m_voxel_size;

            // Count z-layer bounding
            bb.include_z(vs.get_size_z().make(layer_index));

            // Count the bounding
            bb.include_x(vs.get_size_x().make(glom_properties.min_x));
            bb.include_x(vs.get_size_x().make(glom_properties.max_x));
            bb.include_y(vs.get_size_y().make(glom_properties.min_y));
            bb.include_y(vs.get_size_y().make(glom_properties.max_y));
        }
    }

    // Calculate the properties of the glomeruli

    const double glomerulus_min_volume = 4.0 / 3.0 * M_PI * std::pow(m_glomeruli_min_rad.query(), 3);
    const double glomerulus_max_volume = 4.0 / 3.0 * M_PI * std::pow(m_glomeruli_max_rad.query(), 3);
    const double glomerulus_min_diameter = 2 * m_glomeruli_min_rad.query();
    const double glomerulus_max_diameter = 2 * m_glomeruli_max_rad.query();

    for(auto &kv : result.data) {

        glomerulus &glom = kv.second; // TODO: Set location of glomerulus

        glom.label = kv.first;
        glom.volume = glom.pixels.get_volume(module->m_voxel_size);
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

    module->m_output_quantification.attach_foreign(std::move(result), module->m_output_segmented3d);
}

void quantification_constrained_klingberg::create_parameters(misa_parameter_builder &t_parameters) {
    quantification_base::create_parameters(t_parameters);
    m_glomeruli_min_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-min-rad", 15);
    m_glomeruli_max_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-max-rad", 65);
}
