//
// Created by rgerst on 14.12.18.
//

#include <misaxx_kidney_glomeruli/algorithms/quantification/quantification_constrained_klingberg.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;
using namespace coixx;

void quantification_constrained_klingberg::misa_work() {

    for(size_t layer_index = 0; layer_index < m_input_segmented3d.size(); ++layer_index) {

        images::grayscale32s img_components = labels_file->load();

        using namespace coixx::toolbox;
        objects::label_properties<label_pixel_count, label_min_max_position> component_properties(img_components);

        for(const auto& [group, glom_properties] : component_properties) {

            if(group == 0)
                continue;

            auto &glom = module().get_glomerulus(group);
            glom.pixels += glom_properties.get<label_pixel_count>().pixels;

            auto &bb = glom.bounds;
            const auto &lbl_minmax = glom_properties.get<label_min_max_position>();

            // Count z-layer bounding
            bb.min_z = std::min<int>({ bb.min_z, static_cast<int>(layer_index) });
            bb.max_z = std::max<int>({ bb.max_z, static_cast<int>(layer_index) });

            // Count the bounding
            bb.min_x = std::min<int>({ bb.min_x, lbl_minmax.min_x });
            bb.max_x = std::max<int>({ bb.max_x, lbl_minmax.max_x });
            bb.min_y = std::min<int>({ bb.min_y, lbl_minmax.min_y });
            bb.max_y = std::max<int>({ bb.max_y, lbl_minmax.max_y });

        }
    }

    // Calculate the properties of the glomeruli

    const double glomerulus_min_volume = 4.0 / 3.0 * M_PI * pow(m_glomeruli_min_rad, 3);
    const double glomerulus_max_volume = 4.0 / 3.0 * M_PI * pow(m_glomeruli_max_rad, 3);
    const double glomerulus_min_diameter = 2 * m_glomeruli_min_rad;
    const double glomerulus_max_diameter = 2 * m_glomeruli_max_rad;

    for(auto &kv : module().get_glomeruli().data) {

        glomerulus &glom = kv.second;

        glom.label = kv.first;
        glom.volume = glom.pixels * m_voxel_size.volume();
        glom.diameter = 2 * pow(3.0 / 4.0 * glom.volume / M_PI, 1.0 / 3.0);

        const auto &bb = glom.bounds;

        if(bb.is_valid()) {
            const auto sz = bb.to_size();
            const double feret_x = sz.x * m_voxel_size.x;
            const double feret_y = sz.y * m_voxel_size.z;
            const double feret_z = sz.z * m_voxel_size.y;

            glom.valid = glom.volume >= glomerulus_min_volume && glom.volume <= glomerulus_max_volume &&
                         bb.is_valid() &&
                         feret_x >= glomerulus_min_diameter && feret_x <= glomerulus_max_diameter &&
                         feret_y >= glomerulus_min_diameter && feret_y <= glomerulus_max_diameter &&
                         feret_z >= glomerulus_min_diameter && feret_z <= glomerulus_max_diameter;
        }
        else {
            glom.valid = false;
        }
    }
}
