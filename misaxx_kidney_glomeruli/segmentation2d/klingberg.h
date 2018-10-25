//
// Created by rgerst on 11.09.18.
//


#pragma once

#include "segmentation2d_base.h"
#include <misaxx/vdata/object3d_voxel_size.h>

namespace misaxx::module::kidney_glomeruli_detection::segmentation2d {
    struct klingberg : public segmentation2d_base {

        using segmentation2d_base::segmentation2d_base;

        int m_median_filter_size = from_algorithm_json_or<int>("median-filter-size", 3);
        double  m_glomeruli_min_rad = from_algorithm_json_or<double>("glomeruli-min-rad", 15);
        double  m_glomeruli_max_rad = from_algorithm_json_or<double>("glomeruli-max-rad", 65);
        double m_threshold_percentile = from_algorithm_json_or<double>("threshold-percentile", 75);
        double m_threshold_factor = from_algorithm_json_or<double>("threshold-factor", 1.5);
        object3d_voxel_size m_voxel_size = from_parameter<object3d_voxel_size>();

        void misa_work() override {
            using namespace coixx::toolbox;

            images::mask img_non_tissue_mask = m_input_tissue->load();

            if(statistics::is_black(img_non_tissue_mask)) { //INFO: Not inverted yet
                // Instead save a black image
                m_output_segmented2d->save(images::mask(img_non_tissue_mask.get_size(), colors::mask::background()));
                return;
            }

            img_non_tissue_mask << values::invert(); // Now we mark all non-tissue

            images::grayscale_float img = m_input_autofluoresence->auto_load();
            auto img_original = img.clone();

            // Initial median filtering + normalization
            img << blur::median(m_median_filter_size) << normalize::by_max();

            // Generated parameters
            int glomeruli_max_morph_disk_radius = static_cast<int>(m_glomeruli_max_rad / m_voxel_size.xz());
            int glomeruli_min_morph_disk_radius = static_cast<int>((m_glomeruli_min_rad / 2.0) / m_voxel_size.xz());

            // Morphological operation (opening)
            // Corresponds to only allowing objects > disk_size to be included
            // Also subtract the morph result from the initial to remove uneven background + normalize
            img << morph::tophat(structuring_element::ellipse(glomeruli_max_morph_disk_radius * 2)) << normalize::by_max();

            // We are first extracting tissue data
            auto img_tissue = img.clone();

            // Get rid of non-tissue
            img_tissue << values::set_where(colors::grayscale_float::black(), img_non_tissue_mask);

            // Only select glomeruli if the threshold is higher than 75-percentile of kidney tissue
            double percentile_tissue = statistics::get_percentile_as<double>(img, m_threshold_percentile);

            //////////////
            // Now working in uint8
            //////////////

            // Threshold the main image
            uchar otsu_threshold = 0;
            images::grayscale8u img_as8u = semantic_convert<images::mask>(img) << binarize::otsu(otsu_threshold);
            if((otsu_threshold / 255.0) > percentile_tissue * m_threshold_factor ) {

                // Get rid of non-tissue
                img_as8u << values::set_where(colors::grayscale8u::black(), img_non_tissue_mask);

                // Morphological operation (object should have min. radius)
                img_as8u << morph::open(structuring_element::ellipse(glomeruli_min_morph_disk_radius * 2));
            }
            else {
                img_as8u << values::set(colors::grayscale8u(0));
            }

            // Save the mask
            m_output_segmented2d->save(img_as8u);
        }
    };
}