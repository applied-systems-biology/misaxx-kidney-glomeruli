//
// Created by rgerst on 17.12.18.
//

#include "segmentation3d_klingberg.h"
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <utility>
#include <unordered_set>
#include <cv-toolbox/ReadableBMatTypes.h>

using namespace misaxx;
using namespace coixx;
using namespace misaxx_kidney_glomeruli;

/**
* Internally used for the layer graph
*/
struct node_weight {
    /**
     * Layer where the object is located
     */
    size_t layer = 0;

    /**
     * Label within its layer
     */
    int label = 0;

    /**
     * Height of the current object
     */
    size_t height = 0;
};

struct cc_properties {
    void update(int x, int y, int label) {

    }
};

namespace {
    /**
 * Recalculates the height property of each node
 * @param layer_graph
 */
    void recalculate_heights(lemon::ListGraph &g, lemon::ListGraph::NodeMap<node_weight> &w, size_t max_layer) {
        for(auto nd = lemon::ListGraph::NodeIt(g); nd != lemon::INVALID; ++nd) {
            w[nd].height = 0;
        }
        for(auto e = lemon::ListGraph::EdgeIt(g); e != lemon::INVALID; ++e) {
            auto u = g.u(e);
            auto v = g.v(e);
            assert(w[v].layer < w[u].layer);

            w[u].height = max_layer - std::max(w[u].height, w[v].height + 1);
        }
    }

/**
 * Cut all connections to the cutoff-layer from below
 * @param g
 * @param layer_nodes
 * @param cutoff_layer
 */
    void cut_connections_to(lemon::ListGraph &g, const lemon::ListGraph::Node &nd, lemon::ListGraph::NodeMap<node_weight> &w) {
        std::vector<lemon::ListGraph::Edge> to_remove;
        for(auto e = lemon::ListGraph::IncEdgeIt(g, nd); e != lemon::INVALID; ++e) {
            auto u = g.u(e);
            auto v = g.v(e);
            assert(v == nd);

            if(w[u].layer < w[v].layer) {
                to_remove.push_back(e);
            }
        }
        for(const auto &e : to_remove) {
            g.erase(e);
        }
    }
}

void segmentation3d_klingberg::work() {

    // The limit on how large a glomerulus in Z-direction can be at most (this is actually wrong. should be converted from µm to pixels!)
    const size_t maximum_layer_count = static_cast<size_t>(m_max_glomerulus_radius.query());

    // Layers and their names, as well as the number of already saved layers
    cv::images::grayscale32s layer_last;

    // Graph where a node consists of (layer_index, label) and edges represent that
    // two nodes should be assigned to the same final label
    lemon::ListGraph layer_graph;
    lemon::ListGraph::NodeMap<node_weight> node_weights(layer_graph);

    // Assigns the group of the last layer to its LEMON node
    std::vector<std::unordered_map<int, lemon::ListGraph::Node>> layer_nodes;

    // For the first layer, only record the nodes
    {
        layer_nodes.emplace_back(std::unordered_map<int, lemon::ListGraph::Node>());

        const auto input_plane = m_input_segmented2d.at(0);
        auto output_plane = m_output_segmented3d.at(0);

        // Label the 2D segmented object masks
        int img_labels_max_component = 0;
        cv::connectedComponents(input_plane.access_readonly().get(), layer_last, 8, CV_32S);
        output_plane.write(layer_last.clone());

        // Process the components
        cv::label_properties<cc_properties> prop(layer_last);
        for(const auto &kv : prop) {
            auto nd = layer_graph.addNode();
            node_weight o;
            o.layer = 0;
            o.label = kv.first;
            node_weights.set(nd, o);
            layer_nodes[0][kv.first] = nd;
        }

        std::cout << "Found " << (img_labels_max_component - 1) << " glomeruli in first layer" << std::endl;
    }

    // For all other layers also look at the overlap
    for(size_t layer_index = 1; layer_index < m_input_segmented2d.size(); ++layer_index) {
        layer_nodes.emplace_back(std::unordered_map<int, lemon::ListGraph::Node>());

        // Label the 2D segmented object masks
        const auto input_plane = m_input_segmented2d.at(layer_index);
        auto output_plane = m_output_segmented3d.at(layer_index);

        int img_labels_max_component = 0;
        images::grayscale32s img_labels = labeling::connected_components(input_plane.access_readonly().get(), img_labels_max_component);
        output_plane.write(img_labels.clone());

        std::cout << "Found " << img_labels_max_component << " glomeruli in this layer" << std::endl;

        // Find the edges
        std::unordered_set<int> new_nodes;
        std::unordered_set<std::pair<int, int>, boost::hash<std::pair<int, int>>> edges;

        for(int i = 0; i < layer_last.get_mat().rows; ++i) {

            const colors::labels *row_layer_last = layer_last.row_ptr(i);
            const colors::labels *row_layer_current = img_labels.row_ptr(i);

            for(int j = 0; j < layer_last.get_mat().cols; ++j) {

                const int group = row_layer_current[j]; // The object in the current layer
                const int bottom_group = row_layer_last[j]; // The object in the last layer

                if(group > 0) {
                    new_nodes.insert(group);
                }
                if(group > 0 && bottom_group > 0) {
                    edges.insert(std::make_pair(group, bottom_group));
                }
            }
        }

        // Add new nodes into the graph
        for(int u : new_nodes) {
            auto nd = layer_graph.addNode();
            node_weight o;
            o.layer = layer_index;
            o.label = u;
            node_weights.set(nd, o);
            layer_nodes[layer_index][u] = nd;
        }

        // Connect edges
        for(const std::pair<int, int> &uv : edges) {
            int u = uv.first;
            int v = uv.second;
            layer_graph.addEdge(layer_nodes[layer_index].at(u), layer_nodes[layer_index - 1].at(v));
        }

        // Store the last layer from the current groups
        layer_last = img_labels.clone();

//        std::cout << "Layer finished. Current number of non-unique groups is " << layers_group_number  << std::endl;
    }

    // Objects can be still connected with small 1px connections or overdetection
    // Use the height of each of those connected objects to split edges according to a maximum height
    // TODO: In Anna's algorithm maximum_layer_count is equal to int(max glomerulus size), but it should be int(max_glomerulus_size / voxel_size.z)
    // Emulates implementation by Klingberg et al where object updates go from "top to bottom"
    recalculate_heights(layer_graph, node_weights, m_input_segmented2d.size() - 1);
    for(size_t next_layer_index = m_input_segmented2d.size(); next_layer_index != 0; --next_layer_index) {
        for(const auto &kv : layer_nodes.at(next_layer_index - 1)) {
            const auto &nd = kv.second;
            size_t height = (m_input_segmented2d.size() - 1) - node_weights[nd].height; // The height from "top to bottom"
            if(height > maximum_layer_count - 1) {
                cut_connections_to(layer_graph, nd, node_weights);
                recalculate_heights(layer_graph, node_weights, m_input_segmented2d.size() - 1);
            }
        }
    }

    // Find the connected components in the graph and generate a LUT for each layer based on this component
    // Then recolor the layers
    lemon::ListGraph::NodeMap<int> connected_components(layer_graph);
    lemon::connectedComponents(layer_graph, connected_components);
    for(size_t layer_index = 0; layer_index < m_output_segmented3d.size(); ++layer_index) {
        std::cout << "Applying LUT " << std::to_string(layer_index + 1) << " / " << std::to_string(m_output_segmented3d.size()) << std::endl;

        zero_recoloring_hashmap<colors::labels> recoloring;
        for(auto nd = lemon::ListGraph::NodeIt(layer_graph); nd != lemon::INVALID; ++nd) {
            const auto obj = node_weights[nd];
            if(obj.layer == layer_index) {
                int component = connected_components[nd] + 1; // Starts with 0
                recoloring.set_recolor(colors::labels(obj.label), colors::labels(component));
            }
        }

        auto output_plane = m_output_segmented3d.at(layer_index);
        auto rw = output_plane.access_readwrite();
        rw.get() << recoloring::recolor(recoloring);
    }
}

void segmentation3d_klingberg::create_parameters(misa_parameter_builder &t_parameters) {
    segmentation3d_base::create_parameters(t_parameters);
    m_max_glomerulus_radius = t_parameters.create_algorithm_parameter<double>("max-glomerulus-radius", 65);
}
