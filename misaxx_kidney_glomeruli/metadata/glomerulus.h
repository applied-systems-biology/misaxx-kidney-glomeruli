//
// Created by rgerst on 17.09.18.
//


#pragma once

#include <misaxx/parameters/voxel_size.h>
#include <misaxx/metadata/misa_metadata.h>

namespace misaxx::module::kidney_glomeruli_detection {
    struct glomerulus : public misa_metadata {
        /**
         * Number of pixels
         */
        int pixels;
        /**
         * Volume of the glomerulus
         */
        double volume;
        /**
         * Diameter of the glomerulus
         */
        double diameter;
        /**
         * Bounding box of the glomerulus
         */
        voxel_size bounds;
        /**
         * Label in the labeling output
         */
        int label;
        /**
         * True if the glomerulus is detected as valid during quantification
         */
        bool valid;

        nlohmann::json to_json() const override {
            nlohmann::json j;
            j["pixels"] = pixels;
            j["volume"] = volume;
            j["diameter"] = diameter;
            j["bounds"] = bounds;
            j["label"] = label;
            j["valid"] = valid;
            return j;
        }

        std::string get_name() const override {
            return "glomerulus";
        }
    };
}