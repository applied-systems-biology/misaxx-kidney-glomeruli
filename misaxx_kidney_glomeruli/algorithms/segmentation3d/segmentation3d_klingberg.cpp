//
// Created by rgerst on 17.12.18.
//

#include <coixx/toolbox/toolbox_labeling.h>
#include <coixx/recoloring_map.h>
#include <coixx/toolbox/toolbox_recoloring.h>
#include <cxxh/unordered_set_components.h>
#include "segmentation3d_klingberg.h"
//#include <lemon/list_graph.h>

using namespace misaxx;
using namespace coixx;
using namespace misaxx_kidney_glomeruli;

void segmentation3d_klingberg::misa_work() {
    using namespace coixx::toolbox;
    using recoloring_t = identity_recoloring_hashmap<colors::labels>;

    // Generated parameters
    int maximum_layer_count = static_cast<int>(m_max_glomerulus_radius);

    // Tracks the number of groups we know
    int layers_group_number = 0;

    // Layers and their names, as well as the number of already saved layers
    images::grayscale32s layer_last;

    // The LUTs contain assign the glomerulus to the group index
    std::vector<recoloring_t> layer_groups;
    layer_groups.reserve(m_input_segmented2d.size());

    for(size_t layer_index = 0; layer_index < m_input_segmented2d.size(); ++layer_index) {

        const auto input_plane = m_input_segmented2d.at(index);
        auto output_plane = m_output_segmented3d.at(index);

        // Separate the glomeruli objects
        int img_labels_max_component = 0;
        images::grayscale32s img_labels = labeling::connected_components(input_plane.access_readonly().get(), img_labels_max_component);
        output_plane.write(img_labels.clone());

        std::cout << "Found " << (img_labels_max_component - 1) << " glomeruli in this layer" << std::endl;

        if(layer_index > 0) {

            // We iterate through all pixels of last_layer and img_labels and track the groups we encounter
            // The LUT of the previous labels is used to find the glomeruli groups
            std::unordered_map<int, std::unordered_set<int>> encountered_glomeruli_groups;
            std::unordered_set<int> groups_new;

            for(int i = 0; i < layer_last.get_image().rows; ++i) {

                const colors::labels *row_layer_last = layer_last.row_ptr(i);
                const colors::labels *row_layer_current = img_labels.row_ptr(i);

                for(int j = 0; j < layer_last.get_image().cols; ++j) {

                    const int group = row_layer_current[j]; // The group in the current layer
                    const int last_group = row_layer_last[j]; // The RAW (not glomeruli!!!) group in the last layer

                    if(group > 0) {

                        auto it = encountered_glomeruli_groups.find(group);

                        // 2 cases:
                        // 1. Existing group in previous layer: add to set of encountered glomeruli groups
                        // 2. Previous layer has no group and we never encountered the group: Add new group
                        if(last_group > 0) {

                            const int glomeruli_group = layer_groups[layer_groups.size() - 1].recolor(last_group);

                            if(it != encountered_glomeruli_groups.end()) {
                                it->second.insert(glomeruli_group);
                            }
                            else {
                                encountered_glomeruli_groups[group].insert(glomeruli_group);
                            }

                            groups_new.erase(group);
                        }
                        else if(it == encountered_glomeruli_groups.end()){
                            groups_new.insert(group);
                        }
                    }
                }
            }

            // Create a LUT and update it accordingly
            recoloring_t lut;
            lut.data.reserve(groups_new.size() + encountered_glomeruli_groups.size());
//                    lut.reserve(layers_group_number + static_cast<int>(groups_new.size()));

            // Add completely new groups into the LUT
            for(const int group : groups_new) {

                assert(encountered_glomeruli_groups.find(group) == encountered_glomeruli_groups.end());

                if(layers_group_number == std::numeric_limits<int>::max()) {
                    throw std::runtime_error("Reached the maximum number of supported groups!");
                }

                ++layers_group_number;
                lut.set_recolor(group, layers_group_number);
            }

            // Determine the encountered group connected components
            // This will be later used to correctly re-assign the groups
            cxxh::unordered_set_components<int> cc;
            for (const auto &encountered : encountered_glomeruli_groups) {
                const int representative = *(encountered.second.begin());

                // Connect encountered groups to the representative, so all of them are connected
                // When we connect to this group later, we obtain a new representative
                for(const auto v : encountered.second) {
                    cc.insert(v);
                    cc.connect(v, representative);
                }
            }

            // Assign a new group to all already "covered" groups. Needed to allow vertical restriction of labeling
            const auto cc_map = cc.to_map();
            std::unordered_set<int> cc_new;

            for(const auto &kv : cc_map) {
                // Create new group and insert it into the CC
                ++layers_group_number;
                cc.insert(layers_group_number);
                cc.connect(*(kv.second.begin()), layers_group_number);
                cc_new.insert(layers_group_number);
            }

            // We build a lut that transforms the existing groups into the new groups
            const auto cc_lut = cc.to_sink_transform_function<recoloring_t>(cc_new);

            // Re-assign the local groups using the sink_lut
            for (const auto &encountered : encountered_glomeruli_groups) {
                const int raw_group = encountered.first;
                const int target_global_group = cc_lut.recolor(*(encountered.second.begin())); // Important: Use the LUT for joined grouping!

                lut.set_recolor(raw_group, target_global_group);
            }

            for(int i = std::max(0, static_cast<int>(layer_index) - maximum_layer_count); i < static_cast<int>(layer_groups.size()); ++i) {
                auto &previous_lut = layer_groups[i];
                cxxh::math::function::filter(previous_lut, cxxh::math::function::chain_discrete_with(cc_lut));
            }

            layer_groups.push_back(lut);

        }
        else {
            // The first layer is trivial
            // Trivial LUT and unchanged group image
            recoloring_t lut;
            lut.data.reserve(img_labels_max_component);
            layer_groups.push_back(lut);

            layers_group_number = img_labels_max_component;
        }

        // Store the last layer from the current groups
        layer_last = img_labels.clone();

        std::cout << "Layer finished. Current number of non-unique groups is " << layers_group_number  << std::endl;
    }

    // We recolor all labels according to the LUT
    for(size_t layer_index = 0; layer_index < m_output_segmented3d.size(); ++layer_index) {

        std::cout << "Applying LUT " << std::to_string(layer_index + 1) << " / " << std::to_string(m_output_segmented3d.size()) << std::endl;

        auto output_plane = m_output_segmented3d.at(layer_index);
        auto rw = output_plane.access_readwrite();
        const auto &lut = layer_groups[layer_index];
        rw.get() << recoloring::recolor(lut);
    }
}
