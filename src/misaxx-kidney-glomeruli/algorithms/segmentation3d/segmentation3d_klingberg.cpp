//
// Created by rgerst on 08.04.19.
//

#include "segmentation3d_klingberg.h"

namespace cv::images {
    using mask = cv::Mat_<uchar>;
    using labels = cv::Mat_<int>;
}

void misaxx_kidney_glomeruli::segmentation3d_klingberg::work() {

    auto module = get_module_as<module_interface>();

    std::vector<cv::images::labels> labels;
    size_t loaded_label_count = 0;
    size_t first_loaded_label_index = 0;

    const auto limsize = static_cast<size_t>(m_max_glomerulus_radius.query());

    int global_max_label = 0;

    for(size_t i = 0; i < module->m_output_segmented2d.size(); ++i) {
        const auto mask_access = module->m_output_segmented2d.at(i).access_readonly();
        cv::images::labels label;
        int max_label = cv::connectedComponents(mask_access.get(), label, 4, CV_32S);

        std::cout << "Found " << max_label << " glomeruli in layer "<< std::to_string(i) << "\n";

        if(labels.size() > 0) {
            std::unordered_map<int, std::unordered_set<int>> connections;

            // Look for connections to the last layer if available
            {
                const cv::images::labels &last_label = labels[i - 1];
                for(int y = 0; y < label.rows; ++y) {
                    const int *row = label[y];
                    const int *last_row = last_label[y];
                    for(int x = 0; x < label.cols; ++x) {
                        if(row[x] > 0) {
                            if(last_row[x] > 0) {
                                connections[row[x]].insert(last_row[x]);
                            }
                            else {
                                connections[row[x]]; // Declare existance
                            }
                        }
                    }
                }
            }


            // Apply re-naming of previous layers (if multiple have been found)
            {
                // Obtain a renaming map
                std::unordered_map<int, int> renaming_map;
                for(auto &kv : connections) {
                    if(kv.second.size() > 1) {
                        int target = *kv.second.begin();
                        for(int source : kv.second) {
                            renaming_map[source] = target;
                        }

                        // Rename for the current layer connection
                        kv.second = { target };
                    }
                }

                // Go through last layers and rename
                for(size_t j = 0; j < labels.size() - 1; ++j) {
                    cv::images::labels &last_label = labels[j];
                    for(int y = 0; y < last_label.rows; ++y) {
                        int *row = last_label[y];
                        for(int x = 0; x < label.cols; ++x) {
                            if(row[x] > 0) {
                                auto it = renaming_map.find(row[x]);
                                if(it != renaming_map.end()) {
                                    row[x] = it->second;
                                }
                            }
                        }
                    }
                }
            }

            // Go through the current layer and rename to target / create new groups
            {
                std::unordered_map<int, int> renaming_map;
                for(auto &kv : connections) {
                    if (kv.second.size() == 1) {
                        renaming_map[kv.first] = *kv.second.begin();
                    }
                    else if(kv.second.empty()) {
                        ++global_max_label;
                        renaming_map[kv.first] = global_max_label;
                    }
                    else {
                        throw std::runtime_error("Incomplete layer assignment!");
                    }
                }

                for(int y = 0; y < label.rows; ++y) {
                    int *row = label[y];
                    for(int x = 0; x < label.cols; ++x) {
                        if(row[x] > 0) {
                            row[x] = renaming_map.at(row[x]);
                        }
                    }
                }
            }

        }
        else {
            // Set global max label to current glomeruli count
            global_max_label = max_label;
        }

        labels.emplace_back(std::move(label));
        ++loaded_label_count;

        while(loaded_label_count > limsize) {
            module->m_output_segmented3d.at(first_loaded_label_index).access_write().set(std::move(labels.at(first_loaded_label_index)));
            --loaded_label_count;
            ++first_loaded_label_index;
        }
    }

    // Move all other labels into the cache
    for(size_t j = first_loaded_label_index; j < labels.size(); ++j) {
        if(labels.at(j).empty())
            throw std::logic_error("Unexpected unloaded label!");
        module->m_output_segmented3d.at(j).access_write().set(std::move(labels.at(j)));
    }
}

void misaxx_kidney_glomeruli::segmentation3d_klingberg::create_parameters(
        misaxx::misa_task::parameter_list &t_parameters) {
    misa_task::create_parameters(t_parameters);
    m_max_glomerulus_radius = t_parameters.create_algorithm_parameter<double>("max-glomerulus-radius", 65);
}
