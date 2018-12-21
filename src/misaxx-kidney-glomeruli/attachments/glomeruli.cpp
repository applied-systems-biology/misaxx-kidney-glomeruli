//
// Created by rgerst on 17.12.18.
//

#include <misaxx-kidney-glomeruli/attachments/glomeruli.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

void glomeruli::from_json(const nlohmann::json &t_json) {
    misaxx::misa_locatable<misaxx::misa_location>::from_json(t_json);
    data = t_json["data"].get<std::unordered_map<int, glomerulus>>();
}

void glomeruli::to_json(nlohmann::json &t_json) const {
    misaxx::misa_locatable<misaxx::misa_location>::to_json(t_json);
    for(const auto &kv : data) {
        kv.second.to_json(t_json["data"][cxxh::to_string(kv.first)]);
    }
}

void glomeruli::to_json_schema(const misaxx::misa_json_schema &t_schema) const {
    misaxx::misa_locatable<misaxx::misa_location>::to_json_schema(t_schema);
    t_schema.resolve("data").declare_required<std::unordered_map<int, glomerulus>>();
}

void glomeruli::build_serialization_id_hierarchy(std::vector<misaxx::misa_serialization_id> &result) const {
    misaxx::misa_locatable<misaxx::misa_location>::build_serialization_id_hierarchy(result);
    result.emplace_back(misaxx::misa_serialization_id("misa_kidney_glomeruli", "attachments/glomeruli"));
}
