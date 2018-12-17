//
// Created by rgerst on 17.09.18.
//


#pragma once

#include <misaxx/misa_serializeable.h>
#include <misaxx_ome/attachments/misa_ome_voxel.h>
#include <misaxx_ome/attachments/misa_ome_pixel_count.h>

namespace misaxx_kidney_glomeruli {

    struct glomerulus : public misaxx::misa_serializeable {
        /**
         * Number of pixels
         */
        misaxx_ome::misa_ome_pixel_count pixels;
        /**
         * Volume of the glomerulus
         */
        misaxx::misa_quantity<double, misaxx_ome::misa_ome_unit_length<3>> volume;
        /**
         * Diameter of the glomerulus
         */
        misaxx::misa_quantity<double, misaxx_ome::misa_ome_unit_length<1>> diameter;
        /**
         * Bounding box of the glomerulus
         */
        misaxx_ome::misa_ome_voxel bounds;
        /**
         * Label in the labeling output
         */
        int label = 0;
        /**
         * True if the glomerulus is detected as valid during quantification
         */
        bool valid = false;

        glomerulus() = default;

        void from_json(const nlohmann::json &j) override;

        void to_json(nlohmann::json &j) const override;

        void to_json_schema(const misaxx::misa_json_schema &t_schema) const override;

    protected:

        void build_serialization_id_hierarchy(std::vector<misaxx::misa_serialization_id> &result) const override;
    };

    inline void to_json(nlohmann::json &j, const glomerulus &p) {
        p.to_json(j);
    }

    inline void from_json(const nlohmann::json &j, glomerulus &p) {
        p.from_json(j);
    }
}