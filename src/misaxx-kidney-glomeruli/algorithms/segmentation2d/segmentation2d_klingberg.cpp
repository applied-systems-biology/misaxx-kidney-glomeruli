//
// Created by rgerst on 17.12.18.
//

#include "segmentation2d_klingberg.h"
#include <cv-toolbox/ReadableBMatTypes.h>
#include <cv-toolbox/toolbox/toolbox_values.h>
#include <cv-toolbox/toolbox/toolbox_semantic_convert.h>
#include <cv-toolbox/toolbox/toolbox_normalize.h>
#include <cv-toolbox/toolbox/toolbox_morph.h>
#include <cv-toolbox/toolbox/toolbox_statistics.h>
#include <cv-toolbox/toolbox/toolbox_binarize.h>
#include <cv-toolbox/structuring_element.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

void segmentation2d_klingberg::work() {

    auto pam = is_parallelizeable_parameter;

    auto module = get_module_as<module_interface>();

    cv::images::mask img_non_tissue_mask { m_input_tissue.clone() };

    if(cv::countNonZero(img_non_tissue_mask) == 0) { //INFO: Not inverted yet
        // Instead save a black image
        m_output_segmented2d.write(cv::images::mask(img_non_tissue_mask.size(), 0));
        return;
    }

    cv::toolbox::invert(img_non_tissue_mask); // Now we mark all non-tissue

    cv::images::grayscale32f img = cv::toolbox::semantic_convert::to_grayscale32f(m_input_autofluoresence.clone());
    auto img_original = img.clone();

    // Initial median filtering + normalization
    cv::medianBlur(img, img.buffer(), m_median_filter_size.query()); img.swap();
    cv::toolbox::normalize::by_max(img);

    // Generated parameters
    const double voxel_xy = module->m_voxel_size.get_size_xy().get_value();
    int glomeruli_max_morph_disk_radius = static_cast<int>(m_glomeruli_max_rad.query() / voxel_xy);
    int glomeruli_min_morph_disk_radius = static_cast<int>((m_glomeruli_min_rad.query() / 2.0) / voxel_xy);

    // Morphological operation (opening)
    // Corresponds to only allowing objects > disk_size to be included
    // Also subtract the morph result from the initial to remove uneven background + normalize
    cv::toolbox::morph::tophat(img, cv::structuring_element::ellipse(glomeruli_max_morph_disk_radius * 2));
    cv::toolbox::normalize::by_max(img);

    // We are first extracting tissue data
    auto img_tissue = img.clone();

    // Get rid of non-tissue
    img_tissue(img_non_tissue_mask) = 0;

    // Only select glomeruli if the threshold is higher than 75-percentile of kidney tissue
    double percentile_tissue = cv::toolbox::statistics::get_percentile(img, m_threshold_percentile.query());

    //////////////
    // Now working in uint8
    //////////////

    // Threshold the main image
    cv::images::grayscale8u img_as8u = cv::toolbox::semantic_convert::to_grayscale8u(img);
    uchar otsu_threshold = cv::toolbox::otsu(img_as8u);
    if((otsu_threshold / 255.0) > percentile_tissue * m_threshold_factor.query() ) {

        // Get rid of non-tissue
        img_as8u(img_non_tissue_mask) = 0;

        // Morphological operation (object should have min. radius)
        cv::toolbox::morph::open(img_as8u, cv::structuring_element::ellipse(glomeruli_min_morph_disk_radius * 2));
    }
    else {
        img_as8u.self() = 0;
    }

    // Save the mask
    m_output_segmented2d.write(std::move(img_as8u));
}

void segmentation2d_klingberg::create_parameters(misa_parameter_builder &t_parameters) {
    segmentation2d_base::create_parameters(t_parameters);
    m_median_filter_size = t_parameters.create_algorithm_parameter<int>("median-filter-size", 3);
    m_glomeruli_min_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-min-rad", 15);
    m_glomeruli_max_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-max-rad", 65);
    m_threshold_percentile = t_parameters.create_algorithm_parameter<double>("threshold-percentile", 75);
    m_threshold_factor = t_parameters.create_algorithm_parameter<double>("threshold-factor", 1.5);
}
