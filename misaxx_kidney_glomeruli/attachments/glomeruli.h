//
// Created by rgerst on 17.09.18.
//


#pragma once

#include <cxxh/string.h>
#include "glomerulus.h"

namespace misaxx_kidney_glomeruli {

    struct glomeruli : public misaxx::misa_serializeable {
        std::unordered_map<int, glomerulus> data;

        void from_json(const nlohmann::json &t_json) override {

        }

        void to_json(nlohmann::json &t_json) const override {
            misa_serializeable::to_json(t_json);
            for(const auto &kv : data) {
                kv.second.to_json(t_json["data"][cxxh::to_string(kv.first)]);
            }
        }

        void to_json_schema(const misaxx::misa_json_schema &t_schema) const override {
            t_schema.resolve("data").declare_required<std::unordered_map<int, glomerulus>>();
        }

    protected:
        void build_serialization_id_hierarchy(std::vector<misaxx::misa_serialization_id> &result) const override {
            misa_serializeable::build_serialization_id_hierarchy(result);
            result.emplace_back(misaxx::misa_serialization_id("misa_kidney_glomeruli", "attachments/glomeruli"));
        }
    };

    inline void to_json(nlohmann::json& j, const glomeruli& p) {
        p.to_json(j);
    }

    inline void from_json(const nlohmann::json& j, glomeruli& p) {
        p.from_json(j);
    }
}