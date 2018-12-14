//
// Created by rgerst on 17.09.18.
//


#pragma once

#include "glomerulus.h"

namespace misaxx_kidney_glomeruli {

    struct glomeruli : public misaxx::misa_serializeable {
        std::unordered_map<int, glomerulus> data;

        nlohmann::json to_json() const override {
            nlohmann::json j;
            for(const auto &kv : data) {
                j[std::to_string(kv.first)] = kv.second.to_json();
            }
            return j;
        }

        std::string get_name() const override {
            return "glomeruli";
        }
    };

    void to_json(nlohmann::json& j, const glomeruli& p) {
        j = p.to_json();
    }

    void from_json(const nlohmann::json& j, glomeruli& p) {
        for(auto it = j.begin(); it != j.end(); ++it) {
            p.data[std::stoi(it.key())] = it.value();
        }
    }
}