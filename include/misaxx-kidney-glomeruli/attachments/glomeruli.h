//
// Created by rgerst on 17.09.18.
//


#pragma once

#include <misaxx/core/attachments/misa_locatable.h>
#include <misaxx/core/attachments/misa_labeled_object_location.h>
#include <misaxx-kidney-glomeruli/attachments/glomerulus.h>

namespace misaxx_kidney_glomeruli {

    struct glomeruli : public misaxx::misa_locatable {
        std::unordered_map<int, glomerulus> data;

        void from_json(const nlohmann::json &t_json) override;

        void to_json(nlohmann::json &t_json) const override;

        void to_json_schema(const misaxx::misa_json_schema &t_schema) const override;

    protected:
        void build_serialization_id_hierarchy(std::vector<misaxx::misa_serialization_id> &result) const override;
    };

    inline void to_json(nlohmann::json& j, const glomeruli& p) {
        p.to_json(j);
    }

    inline void from_json(const nlohmann::json& j, glomeruli& p) {
        p.from_json(j);
    }
}