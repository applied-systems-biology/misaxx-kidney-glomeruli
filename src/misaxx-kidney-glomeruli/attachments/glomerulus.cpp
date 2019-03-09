//
// Created by rgerst on 14.12.18.
//

#include <misaxx-kidney-glomeruli/attachments/glomerulus.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

glomerulus::glomerulus() : volume(misaxx::ome::misa_ome_volume<double> { misaxx::ome::units::micrometer<3>() }),
                           diameter(misaxx::ome::misa_ome_length<double> { misaxx::ome::units::micrometer<1>() }),
                           bounds(misaxx::ome::misa_ome_voxel { misaxx::ome::units::micrometer<1>() }) {
}


void glomerulus::from_json(const nlohmann::json &j) {
    misa_locatable::from_json(j);
    pixels.from_json(j["pixels"]);
    volume.from_json(j["volume"]);
    diameter.from_json(j["diameter"]);
    bounds.from_json(j["bounds"]);
    label = j["label"];
    valid = j["valid"];
}

void glomerulus::to_json(nlohmann::json &j) const {
    misa_locatable::to_json(j);
    pixels.to_json(j["pixels"]);
    volume.to_json(j["volume"]);
    diameter.to_json(j["diameter"]);
    bounds.to_json(j["bounds"]);
    j["label"] = label;
    j["valid"] = valid;
}

void glomerulus::to_json_schema(misaxx::misa_json_schema_property &t_schema) const {
    misa_locatable::to_json_schema(t_schema);
    pixels.to_json_schema(t_schema["pixels"]);
    volume.to_json_schema(t_schema["volume"]);
    diameter.to_json_schema(t_schema["diameter"]);
    bounds.to_json_schema(t_schema["bounds"]);
    t_schema.resolve("label")->declare_required<int>();
    t_schema.resolve("valid")->declare_required<bool>();
}

void glomerulus::build_serialization_id_hierarchy(std::vector<misaxx::misa_serialization_id> &result) const {
    misa_locatable::build_serialization_id_hierarchy(result);
    result.emplace_back(misaxx::misa_serialization_id("misa-kidney-glomeruli", "attachments/glomerulus"));
}

std::string glomerulus::get_documentation_name() const {
    return "Glomerulus";
}

std::string glomerulus::get_documentation_description() const {
    return "Glomerulus detected by the segmentation algorithm";
}

