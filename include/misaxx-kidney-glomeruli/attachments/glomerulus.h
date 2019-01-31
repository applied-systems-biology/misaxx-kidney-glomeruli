//
// Created by rgerst on 17.09.18.
//


#pragma once

#include <misaxx/ome/attachments/misa_ome_voxel.h>
#include <misaxx/ome/attachments/misa_ome_pixel_count.h>
#include <misaxx/ome/attachments/misa_ome_quantity.h>
#include <misaxx/core/attachments/misa_locatable.h>
#include <misaxx/core/attachments/misa_labeled_object_location.h>

namespace misaxx_kidney_glomeruli {

    struct glomerulus : public misaxx::misa_locatable {
        /**
         * Number of pixels
         */
        misaxx::ome::misa_ome_pixel_count pixels;
        /**
         * Volume of the glomerulus
         */
        misaxx::ome::misa_ome_volume<double> volume;
        /**
         * Diameter of the glomerulus
         */
        misaxx::ome::misa_ome_length<double> diameter;
        /**
         * Bounding box of the glomerulus
         */
        misaxx::ome::misa_ome_voxel bounds;
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