//
// Created by rgerst on 17.12.18.
//

#include <misaxx-coixx/toolbox/toolbox_labeling.h>
#include <misaxx-coixx/recoloring_map.h>
#include <misaxx-coixx/toolbox/toolbox_recoloring.h>
#include <misaxx-helpers/unordered_set_components.h>
#include <misaxx-kidney-glomeruli/algorithms/segmentation3d/segmentation3d_klingberg.h>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <misaxx-coixx/toolbox/toolbox_objects.h>
#include <misaxx-coixx/objects/label_dummy_property.h>
#include <utility>

using namespace misaxx;
using namespace coixx;
using namespace misaxx_kidney_glomeruli;

void segmentation3d_klingberg::misa_work() {
    using namespace coixx::toolbox;
    using recoloring_t = identity_recoloring_hashmap<colors::labels>;

    // The limit on how large a glomerulus in Z-direction can be at most
    const auto maximum_layer_count = static_cast<int>(m_max_glomerulus_radius);

    // Layers and their names, as well as the number of already saved layers
    images::grayscale32s layer_last;

    // Graph where a node consists of (layer_index, label) and edges represent that
    // two nodes should be assigned to the same final label
    lemon::ListGraph layer_graph;
    lemon::ListGraph::NodeMap<object> node_map(layer_graph);

    // Assigns the group of the last layer to its LEMON node
    std::unordered_map<int, lemon::ListGraph::Node> last_layer_nodes;
    std::unordered_map<int, lemon::ListGraph::Node> current_layer_nodes;

    // For the first layer, only record the nodes
    {
        const auto input_plane = m_input_segmented2d.at(0);
        auto output_plane = m_output_segmented3d.at(0);

        // Label the 2D segmented object masks
        int img_labels_max_component = 0;
        images::grayscale32s img_labels = labeling::connected_components(input_plane.access_readonly().get(), img_labels_max_component);
        output_plane.write(img_labels.clone());

        // Process the components
        toolbox::objects::label_properties<label_dummy_property> prop(img_labels);
        for(const auto &kv : prop) {
            auto nd = layer_graph.addNode();
            object o;
            o.layer = 0;
            o.label = kv.first;
            node_map.set(nd, o);
            current_layer_nodes[kv.first] = nd;
        }

        std::cout << "Found " << (img_labels_max_component - 1) << " glomeruli in first layer" << std::endl;
    }

    // For all other layers also look at the overlap
    for(size_t layer_index = 1; layer_index < m_input_segmented2d.size(); ++layer_index) {

        // Clear the node assignments of the current map
        std::swap(current_layer_nodes, last_layer_nodes);
        current_layer_nodes.clear();

        // Label the 2D segmented object masks
        const auto input_plane = m_input_segmented2d.at(layer_index);
        auto output_plane = m_output_segmented3d.at(layer_index);

        int img_labels_max_component = 0;
        images::grayscale32s img_labels = labeling::connected_components(input_plane.access_readonly().get(), img_labels_max_component);
        output_plane.write(img_labels.clone());

        std::cout << "Found " << (img_labels_max_component - 1) << " glomeruli in this layer" << std::endl;

        // Find the edges
        std::unordered_set<int> new_nodes;
        std::unordered_set<std::pair<int, int>, boost::hash<std::pair<int, int>>> edges;

        for(int i = 0; i < layer_last.get_image().rows; ++i) {

            const colors::labels *row_layer_last = layer_last.row_ptr(i);
            const colors::labels *row_layer_current = img_labels.row_ptr(i);

            for(int j = 0; j < layer_last.get_image().cols; ++j) {

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
            object o;
            o.layer = layer_index;
            o.label = u;
            node_map.set(nd, o);
            current_layer_nodes[u] = nd;
        }

        // Connect edges TODO: Anna's algorithm with a limit
        for(const std::pair<int, int> &uv : edges) {
            int u = uv.first;
            int v = uv.second;
            layer_graph.addEdge(current_layer_nodes.at(u), last_layer_nodes.at(v));
        }

        // Store the last layer from the current groups
        layer_last = img_labels.clone();

//        std::cout << "Layer finished. Current number of non-unique groups is " << layers_group_number  << std::endl;
    }

    // Find the connected components in the graph and generate a LUT for each layer based on this component
    // Then recolor the layers
    lemon::ListGraph::NodeMap<int> connected_components(layer_graph);
    lemon::connectedComponents(layer_graph, connected_components);
    for(size_t layer_index = 0; layer_index < m_output_segmented3d.size(); ++layer_index) {
        std::cout << "Applying LUT " << std::to_string(layer_index + 1) << " / " << std::to_string(m_output_segmented3d.size()) << std::endl;

        zero_recoloring_hashmap<colors::labels> recoloring;
        for(auto nd = lemon::ListGraph::NodeIt(layer_graph); nd != lemon::INVALID; ++nd) {
            const auto obj = node_map[nd];
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
