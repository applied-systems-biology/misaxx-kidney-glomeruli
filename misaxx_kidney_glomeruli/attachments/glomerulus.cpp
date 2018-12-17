//
// Created by rgerst on 14.12.18.
//

#include <misaxx_kidney_glomeruli/attachments/glomerulus.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

void glomerulus::from_json(const nlohmann::json &j) {
    pixels.from_json(j["pixels"]);
    volume.from_json(j["volume"]);
    diameter.from_json(j["diameter"]);
    bounds.from_json(j["bounds"]);
    label = j["label"];
    valid = j["valid"];
}

void glomerulus::to_json(nlohmann::json &j) const {
    misa_serializeable::to_json(j);
    pixels.to_json(j["pixels"]);
    volume.to_json(j["volume"]);
    diameter.to_json(j["diameter"]);
    bounds.to_json(j["bounds"]);
    j["label"] = label;
    j["valid"] = valid;
}

void glomerulus::to_json_schema(const misaxx::misa_json_schema &t_schema) const {
    pixels.to_json_schema(t_schema.resolve("pixels"));
    volume.to_json_schema(t_schema.resolve("volume"));
    diameter.to_json_schema(t_schema.resolve("diameter"));
    bounds.to_json_schema(t_schema.resolve("bounds"));
    t_schema.resolve("label").declare_required<int>();
    t_schema.resolve("valid").declare_required<bool>();
}

void glomerulus::build_serialization_id_hierarchy(std::vector<misaxx::misa_serialization_id> &result) const {
    misa_serializeable::build_serialization_id_hierarchy(result);
    result.emplace_back(misaxx::misa_serialization_id("misa_kidney_glomeruli", "attachments/glomerulus"));
}
