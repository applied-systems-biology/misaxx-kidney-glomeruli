//
// Created by rgerst on 08.04.19.
//

#include "segmentation3d_klingberg.h"
#include <set>

namespace cv::images {
    using mask = cv::Mat_<uchar>;
    using labels = cv::Mat_<int>;
}

namespace {
    template<typename Set1, typename Set2>
    bool set_intersects(const Set1 & set1, const Set2 &set2) {
        for(const auto &entry : set1) {
            if(set2.find(entry) != set2.end())
                return true;
        }
        return false;
    }
}

void misaxx_kidney_glomeruli::segmentation3d_klingberg::work() {

    auto module = get_module_as<module_interface>();

    std::vector<cv::images::labels> labels;
    size_t loaded_label_count = 0;
    size_t first_loaded_label_index = 0;

    const auto limsize = static_cast<size_t>(m_max_glomerulus_radius.query());

    int global_max_label = 0;

    for(size_t i = 0; i < module->m_output_segmented2d.size(); ++i) {
        cv::images::labels label;
        int max_label = cv::connectedComponents( module->m_output_segmented2d.at(i).access_readonly().get(),
                label, 4, CV_32S);

        std::cout << "Found " << max_label << " glomeruli in layer "<< std::to_string(i) << "\n";

        if(!labels.empty()) {
            // All connections from this layer -> labels of last layer
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

            // Go through the current layer and rename to target / create new groups
            // The first pass ignores merging objects, but puts the current layer into the
            std::unordered_map<int, int> first_pass_renaming;
            {
                for(auto &kv : connections) {
                    if (kv.second.empty()) {
                        ++global_max_label;
                        first_pass_renaming[kv.first] = global_max_label;
                    }
                    else {
                        first_pass_renaming[kv.first] = *kv.second.begin();
                    }
                }

                for(int y = 0; y < label.rows; ++y) {
                    int *row = label[y];
                    for(int x = 0; x < label.cols; ++x) {
                        if(row[x] > 0) {
                            row[x] = first_pass_renaming.at(row[x]);
                        }
                    }
                }
            }

            // Move the first pass results into the current buffer
            labels.emplace_back(std::move(label));
            ++loaded_label_count;

            // First-pass renaming does not merge objects
            // We need a second renaming pass that merges all connected objects
            {
                // Find the connected components of connections
                std::vector<std::unordered_set<int>> components;
                for(const auto &kv : connections) {
                    if(kv.second.size() > 1)
                        components.push_back(kv.second);
                }

                bool changed;
                do {
                    changed = false;
                    for(size_t j = 0; j < components.size(); ++j) {
                        if(j >= components.size() - 1)
                            continue;
                        for(size_t k = 0; k < components.size(); ++k) {
                            if(j != k && k < components.size()) {
                                if(set_intersects(components.at(j), components.at(k))) {
                                    // Copy the labels from k to j
                                    std::copy(components.at(k).begin(), components.at(k).end(),
                                            std::inserter(components.at(j), components.at(j).begin()));

                                    // Remove k
                                    std::swap(components.at(k), components.at(components.size() - 1));
                                    components.resize(components.size() - 1);


                                    changed = true;
                                }
                            }
                        }
                    }
                }
                while(changed);

                // Create a second pass renaming
                std::unordered_map<int, int> second_pass_renaming;
                for(const auto &component : components) {
                    int target = *component.begin();
                    for(int src : component) {
                        second_pass_renaming[src] = target;
                    }
                }

                // Go through last layers and rename
                // This includes the current layer
                for(size_t j = first_loaded_label_index; j < labels.size(); ++j) {
                    cv::images::labels &current_label = labels[j];

                    for(int y = 0; y < current_label.rows; ++y) {
                        int *row = current_label[y];
                        for(int x = 0; x < current_label.cols; ++x) {
                            if(row[x] > 0) {
                                auto it = second_pass_renaming.find(row[x]);
                                if(it != second_pass_renaming.end()) {
                                    row[x] = it->second;
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            // Set global max label to current glomeruli count
            global_max_label = max_label;
            labels.emplace_back(std::move(label));
            ++loaded_label_count;
        }

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
