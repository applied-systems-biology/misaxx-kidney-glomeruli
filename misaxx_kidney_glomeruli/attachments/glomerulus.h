//
// Created by rgerst on 17.09.18.
//


#pragma once

#include <misaxx/misa_serializeable.h>

namespace misaxx_kidney_glomeruli {

struct glomerulus : public misaxx::misa_serializeable {
        /**
         * Number of pixels
         */
        int pixels = 0;
        /**
         * Volume of the glomerulus
         */
        double volume = 0;
        /**
         * Diameter of the glomerulus
         */
        double diameter = 0;
        /**
         * Bounding box of the glomerulus
         */
        object3d_voxel_bounds bounds;
        /**
         * Label in the labeling output
         */
        int label = 0;
        /**
         * True if the glomerulus is detected as valid during quantification
         */
        bool valid = false;

        nlohmann::json to_json() const override {
            nlohmann::json j;
            j["pixels"] = pixels;
            j["volume"] = volume;
            j["diameter"] = diameter;
            j["bounds"] = bounds.to_json();
            j["label"] = label;
            j["valid"] = valid;
            return j;
        }

        std::string get_name() const override {
            return "glomerulus";
        }


    };

    void to_json(nlohmann::json& j, const glomerulus& p) {
        j = p.to_json();
    }

    void from_json(const nlohmann::json& j, glomerulus& p) {
        p.pixels = j["pixels"];
        p.volume = j["volume"];
        p.diameter = j["diameter"];
        p.bounds = j["bounds"];
        p.label = j["label"];
        p.valid = j["valid"];
    }
}