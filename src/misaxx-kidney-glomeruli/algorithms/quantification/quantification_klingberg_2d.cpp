//
// Created by ruman on 23.03.19.
//

#include "quantification_klingberg_2d.h"
#include "quantification_klingberg.h"
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

        void update(int x, int y, int label) {
            ++pixels;
        }
    };

}


void quantification_klingberg_2d::work() {
    auto module = get_module_as<module_interface>();

    if(m_input_segmented3d.size() != 1)
        throw std::logic_error("Quantification method assumes 2D data!");
    if(module->m_voxel_size.get_size_z().get_value() != 1)
        throw std::logic_error("Voxel depth must be 1!");

    glomeruli result;

    for(const auto &plane : m_input_segmented3d) {
        auto access = plane.access_readonly();

        cv::label_properties<cc_properties> component_properties(cv::images::labels { access.get() });

        for(const auto& [group, glom_properties] : component_properties.rows) {

            if(group == 0)
                continue;

            glomerulus &glom = result.data[group];
            glom.pixels.count += glom_properties.pixels;
        }
    }

    // Calculate the properties of the glomeruli
    double glomerulus_min_area = M_PI * std::pow(m_glomeruli_min_rad.query(), 2);
    double glomerulus_max_area = M_PI * std::pow(m_glomeruli_max_rad.query(), 2);

    double diameter_sum = 0;
    double diameter_sum_sq = 0;

    for(auto &kv : result.data) {
        glomerulus &glom = kv.second;
        glom.label = kv.first;
        glom.volume = glom.pixels.get_volume(module->m_voxel_size); // Equal to area here
        glom.diameter = misa_quantity<double, misa_ome_unit_length<1>>(2 * pow(glom.volume.get_value() / M_PI, 1.0 / 2.0),
                                                                       misa_ome_unit_length<1>::ome_unit_type::MICROMETER);
        glom.valid = glom.volume.get_value() >= glomerulus_min_area && glom.volume.get_value() <= glomerulus_max_area;
        if(glom.valid) {
            ++result.valid_glomeruli_number;
            diameter_sum += glom.diameter.get_value();
            diameter_sum_sq += std::pow(glom.diameter.get_value(), 2);
        }
        else {
            ++result.invalid_glomeruli_number;
        }
    }

    result.valid_glomeruli_diameter_average = diameter_sum / result.valid_glomeruli_number;
    result.valid_glomeruli_diameter_variance = (diameter_sum_sq / result.valid_glomeruli_number) -
                                               std::pow(result.valid_glomeruli_diameter_average, 2);

    module->m_output_quantification.attach_foreign(std::move(result), module->m_output_segmented3d);
}

void quantification_klingberg_2d::create_parameters(misa_parameter_builder &t_parameters) {
    quantification_base::create_parameters(t_parameters);
    m_glomeruli_min_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-min-rad", 15);
    m_glomeruli_max_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-max-rad", 65);
}
