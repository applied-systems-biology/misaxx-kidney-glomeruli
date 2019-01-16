//
// Created by rgerst on 17.12.18.
//

#include "segmentation2d_local_otsu.h"
#include <misaxx-coixx/toolbox/toolbox_normalize.h>
#include <misaxx-coixx/toolbox/toolbox_labeling.h>
#include <misaxx-coixx/toolbox/toolbox_bitwise.h>
#include <misaxx-coixx/toolbox/toolbox_resize.h>
#include <misaxx-coixx/toolbox/toolbox_recoloring.h>
#include <misaxx-coixx/toolbox/toolbox_mask.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;
using namespace coixx;

void segmentation2d_local_otsu::work() {
    using namespace coixx::toolbox;

    auto module = get_module_as<kidney_glomeruli>();

    images::mask tissue_mask = m_input_tissue.clone();

    if(statistics::is_black(tissue_mask)) { //INFO: Not inverted yet
        // Instead save a black image
        images::mask img(tissue_mask.get_size(), colors::mask::background());
        m_output_segmented2d.write(std::move(img));
        return;
    }

    images::grayscale_float img = m_input_autofluoresence.clone();
    images::grayscale_float img_original;

    // Smooth & normalize
    img << values::backup(img_original) << blur::median(m_median_filter_size.query()) << normalize::by_max();

    // Extracts the blobs
    images::grayscale_float blobs = extract_blobs_log(img);

    images::mask cortex_mask = find_cortex_otsu_distance_and_dilation(tissue_mask, blobs);
//            image_filter2::visualization::mask(cortex_mask, img).show_and_wait("cortex visualization");


    img << values::set_where_not(colors::grayscale_float::black(), cortex_mask);

    // Run multiple iterations of local Otsu and merge them if there are no conflicts
    images::mask merged_mask(img.get_size(), colors::mask::background());

    // Local maxima selection is expensive (Due to dilation).
    // Instead use a two-step approach that only requires the dilation with the small selection
    const double voxel_xy = module->m_voxel_size.get_size_xy().get_value();
    images::mask blobs_all_maxima = extract_maxima(blobs.clone(), tissue_mask, m_glomeruli_min_rad.query() / voxel_xy);
    blobs_all_maxima << values::set_where_not(colors::mask::background(), cortex_mask);

    // We first try to segment large glomeruli to prevent oversegmentation (in combination with deleting already segmented maxima)
    for(const double radius_microns : { m_glomeruli_max_rad.query(), m_glomeruli_min_rad.query() }) {
        const double radius = radius_microns / voxel_xy;

        images::mask blobs_maxima = blobs_all_maxima.clone();
        if(radius_microns != m_glomeruli_min_rad.query()) {
            restrict_maxima(blobs_maxima, img, radius);
        }

        images::grayscale32s blobs_maxima_components = labeling::connected_components<colors::grayscale32s >(blobs_maxima);
        images::mask final_mask = segment_glomeruli_local_otsu(blobs_maxima, cortex_mask, img);

        // Delete only good positions from the list of maxima to be analyzed.
        // The reason behind this is that a larger search radius can cover two adjacent glomeruli where
        // only one maximum is seen as relevant. This will lead to bad detection. But if the maximum is already
        // deleted, the small scale algorithm won't be able to detect the right glomeruli.
        blobs_all_maxima << values::set_where(colors::mask::background(), final_mask);

        merged_mask << bitwise::OR(final_mask);
    }

//            image_filter2::visualization::mask(merged_mask, img).show_and_wait("final mask");
    m_output_segmented2d.write(std::move(merged_mask));
}

coixx::images::grayscale_float segmentation2d_local_otsu::extract_blobs_log(const coixx::images::grayscale_float &img) {

    auto module = get_module_as<kidney_glomeruli>();

    using namespace coixx::toolbox;

    const double voxel_xy = module->m_voxel_size.get_size_xy().get_value();
    const double glomeruli_min_rad_sigma = (m_glomeruli_min_rad.query() / voxel_xy) / sqrt(2);
    const double glomeruli_max_rad_sigma = (m_glomeruli_max_rad.query() / voxel_xy) / sqrt(2);
    const double avg_rad_sigma = (glomeruli_min_rad_sigma + glomeruli_max_rad_sigma) / 2;

    images::grayscale_float log_response = img.clone() << blob::laplacian_of_gaussian(avg_rad_sigma);

    // Manual postprocessing for increased speed:
    // Takes only the negative response and normalizes against the lowest response (where the blobs are)
    float min_response = statistics::min(log_response);
    for(int y = 0; y < img.get_image().rows; ++y) {
        colors::grayscale_float *row = log_response.row_ptr(y);
        for(int x = 0; x < img.get_image().cols; ++x) {
            row[x].value = std::min(0.0f, row[x].value) / min_response;
        }
    }

    return log_response;
}

coixx::images::mask segmentation2d_local_otsu::extract_maxima(coixx::images::grayscale_float blobs,
                                                              const coixx::images::mask &tissue_mask,
                                                              const double t_radius) {

    using namespace coixx::toolbox;

    // Find the local maxima
    images::mask blobs_mask = localminmax::local_exclusive_max_morph(blobs, static_cast<int>(2 * t_radius));
    blobs_mask << values::set_where_not(colors::mask::black(), tissue_mask);

    return blobs_mask;
}

void
segmentation2d_local_otsu::restrict_maxima(coixx::images::mask &t_maxima, const coixx::images::grayscale_float &img,
                                           const double t_radius) {

    using namespace coixx::toolbox;

    std::vector<pixel<colors::grayscale_float >> pixels;
    img << values::get_where(std::back_inserter(pixels), t_maxima);

    const double r_sq = pow(2 * t_radius, 2);

    // For each pixel, check if we have another pixel with same or larger value within radius
    // If this is the case, we drop the pixel
    for(size_t i = 0; i < pixels.size(); ++i) {
        const auto px = pixels[i];
        for(size_t j = 0; j < pixels.size(); ++j) {
            if(i != j) {
                const auto px2 = pixels[j];
                const double l = pow(px.x - px2.x, 2) + pow(px.y - px2.y, 2);
                if(l <= r_sq && px2.value >= px.value) {
                    t_maxima.get_image().at<uchar>(cv::Point(px.x, px.y)) = 0; // Get rid of the maximum
                    break;
                }
            }
        }
    }
}

coixx::images::mask
segmentation2d_local_otsu::find_cortex_otsu_distance_and_dilation(const coixx::images::mask &tissue_mask,
                                                                  const coixx::images::grayscale_float &blobs) {

    auto module = get_module_as<kidney_glomeruli>();

    using namespace coixx::toolbox;

    const double voxel_xy = module->m_voxel_size.get_size_xy().get_value();
    const int glomeruli_max_diameter = static_cast<int>((m_glomeruli_max_rad.query() / voxel_xy) * 2);

    images::grayscale_float cortex_dist(blobs.get_size(), colors::grayscale_float::black());
    cv::distanceTransform(tissue_mask.get_image(), cortex_dist.get_image(), cv::DIST_L2, 3);

    // Apply otsu
    images::grayscale8u blobs_thresholded = semantic_convert<images::grayscale8u >(blobs) << binarize::otsu_where(tissue_mask);

    // Dilation based method
    int dilate_diameter = static_cast<int>(glomeruli_max_diameter * 0.1 * m_cortex_segmentation_dilation_group_size.query());
    images::grayscale8u blobs_thresholded_small = resize(blobs_thresholded, 0.1, resize_interpolation::nearest);
    blobs_thresholded_small << morph::dilate(structuring_element::ellipse(dilate_diameter));

    images::mask cortex_by_dilation = resize(blobs_thresholded_small,
                                             blobs_thresholded.get_size(),
                                             resize_interpolation::nearest);
    cortex_by_dilation << values::set_where_not(colors::mask(0), tissue_mask);
    blobs_thresholded << binarize::otsu();

    // Narrow down to glomeruli regions & tissue
    auto cortex_dist_candidates = cortex_dist.clone();
    cortex_dist_candidates << values::set_where_not(colors::grayscale_float(0), tissue_mask);
    cortex_dist_candidates << values::set_where_not(colors::grayscale_float(0), blobs_thresholded);

    // Options: Maximum (works well in test images, BUT outliers disturb it!)
    // Alternative: 2 * mean (seems to work)
//            float dist_max = toolbox(cortex_dist_candidates).max();
    images::mask cortex_dist_candidates_nonzero(cortex_dist_candidates.get_image() > 0);
    double dist = statistics::get_mean_where(cortex_dist_candidates, cortex_dist_candidates_nonzero) * 2;

    images::mask cortex(cortex_dist.get_image() <= dist);

    // Combine both
    cortex << values::set_where(colors::mask(255), cortex_by_dilation);
    cortex << values::set_where_not(colors::mask(0), tissue_mask);

    return cortex;
}

void
segmentation2d_local_otsu::segment_glomeruli_local_otsu_blacklist_by_contour(const coixx::images::grayscale_float &,
                                                                             const coixx::images::mask &blobs_maxima_mask,
                                                                             const coixx::images::mask &local_otsu_mask,
                                                                             const coixx::images::grayscale32s &local_otsu_mask_components,
                                                                             const int local_otsu_mask_max_component_id,
                                                                             coixx::mutable_recoloring_map<coixx::colors::labels> &blacklist) {

    auto module = get_module_as<kidney_glomeruli>();

    using namespace coixx::toolbox;

//            visualizer_gui() << named_image(img, "original") << named_image(blobs_maxima_mask, "blobs maxima") << named_image(local_otsu_mask, "local otsu mask")
//                                                                                                               << named_image(local_otsu_mask_components, "mask components") << show_visualizer_gui();

    // Detect the contours
    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(local_otsu_mask.get_image(), contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    // Find the center for each label
    std::vector<cv::Point> label_centers;
    label_centers.resize(static_cast<size_t>(local_otsu_mask_max_component_id) + 1, cv::Point(-1, -1));

    for(int y = 0; y < local_otsu_mask_components.get_image().rows; ++y) {

        const colors::mask * row_maxima = blobs_maxima_mask.row_ptr(y);
        const colors::labels * row_components = local_otsu_mask_components.row_ptr(y);

        for(int x = 0; x < local_otsu_mask_components.get_image().cols; ++x) {
            if(row_maxima[x] > 0 && row_components[x] > 0) {
                label_centers[row_components[x]] = cv::Point(x, y);
            }
        }
    }

    // Assign a contour to each label
    std::vector<std::vector<cv::Point>> label_contours;
    label_contours.resize(label_centers.size());

    for(auto &contour : contours) {
        for(size_t label = 1; label < label_centers.size(); ++label) {
            if(label_contours[label].empty() && label_centers[label].x >= 0) {
                double d = cv::pointPolygonTest(contour, label_centers[label], false);
                if(d > 0) {
                    label_contours[label] = std::move(contour);
                    break;
                }
            }
        }
    }

    // Blacklist phase
    const double voxel_xy = module->m_voxel_size.get_size_xy().get_value();
    const int glomeruli_min_rad = static_cast<int>(m_glomeruli_min_rad.query() / voxel_xy);
    const int glomeruli_max_rad = static_cast<int>(m_glomeruli_max_rad.query() / voxel_xy);

    for(size_t label = 1; label < label_contours.size(); ++label) {
        const auto contour = label_contours[label];
        const auto center = label_centers[label];

        // We need min 5 points to fit the ellipse
        if(contour.size() >= 5 && center.x >= 0) {

            const auto fitted_ellipse = cv::fitEllipse(contour);
            const double fitted_ellipse_r1 = fitted_ellipse.size.width / 2.0;
            const double fitted_ellipse_r2 = fitted_ellipse.size.height / 2.0;

            // Radius check
            double min_rad = std::min(fitted_ellipse_r1, fitted_ellipse_r2);
            double max_rad = std::max(fitted_ellipse_r1, fitted_ellipse_r2);

            if(std::floor(min_rad) < glomeruli_min_rad ||
               std::floor(max_rad) > glomeruli_max_rad ||
               std::floor(max_rad) < glomeruli_min_rad) {
                blacklist.set_recolor(colors::labels(label), colors::labels::background());
//                        std::cout << "Blacklist " << label << " failed radius requirements" << std::endl;
                continue;
            }

            // Eccentricity check
            if(max_rad / min_rad >= 2) {
                blacklist.set_recolor(colors::labels(label), colors::labels::background());
//                        std::cout << "Blacklist " << label << " failed eccentricity requirements" << std::endl;
                continue;
            }

            // Countour roundness check
            const double isoperimetric_quotient = (4 * M_PI * cv::contourArea(contour)) / pow(cv::arcLength(contour, true), 2);
            const double fitted_ellipse_area = M_PI * fitted_ellipse_r1 * fitted_ellipse_r2;
            const double fitted_ellipse_perimeter = 2 * M_PI * sqrt((pow(fitted_ellipse_r1, 2) + pow(fitted_ellipse_r2, 2)) / 2.0);
            const double fitted_isoperimetric_quotient = (4 * M_PI * fitted_ellipse_area) / pow(fitted_ellipse_perimeter, 2);

//                    std::cout << "Q = " << isoperimetric_quotient << ", expected Q = " << fitted_isoperimetric_quotient << std::endl;

            if(isoperimetric_quotient < m_isoperimetric_quotient_threshold.query() * fitted_isoperimetric_quotient || isoperimetric_quotient > 1) {
                blacklist.set_recolor(colors::labels(label), colors::labels::background());
//                        std::cout << "Blacklist " << label << " failed isoperimetric quotient requirements" << std::endl;
                continue;
            }

        }
        else {
//                    std::cout << "Blacklist " << label << " has no contour" << std::endl;
            blacklist.set_recolor(colors::labels(label), colors::labels::background());
        }
    }

}

coixx::images::mask segmentation2d_local_otsu::segment_glomeruli_local_otsu(const coixx::images::mask &blobs_maxima,
                                                                            const coixx::images::mask &cortex_mask,
                                                                            const coixx::images::grayscale_float &img) {
    auto module = get_module_as<kidney_glomeruli>();

    using namespace coixx::toolbox;

    const double voxel_xy = module->m_voxel_size.get_size_xy().get_value();
    const int glomeruli_max_diameter = static_cast<int>((m_glomeruli_max_rad.query() / voxel_xy) * 2);
    const int glomeruli_search_diameter = glomeruli_max_diameter + m_voronoi_cell_radius_border.query();

    // Create an area around each maximum
    images::mask local_max_areas_mask = blobs_maxima.clone() << morph::dilate(structuring_element::ellipse(glomeruli_search_diameter));
    local_max_areas_mask << values::set_where_not(colors::mask::background(), cortex_mask); // Remove regions outside cortex

    // Generate voronoi partitioning
    images::grayscale_float local_areas_voronoi_distance(img.get_size(), colors::grayscale_float::black());
    images::grayscale8u local_areas_voronoi_centers = blobs_maxima.clone() << values::invert();
    images::grayscale32s local_areas_voronoi_labels(img.get_size(), colors::labels::background());
    cv::distanceTransform(local_areas_voronoi_centers.get_image(),
                          local_areas_voronoi_distance.get_image(),
                          local_areas_voronoi_labels.get_image(),
                          cv::DIST_L2,
                          3);

    // Delete partitioning outside cortex
    local_areas_voronoi_labels << values::set_where_not(colors::grayscale32s::background(), cortex_mask);

    // Apply per-component localized Otsu on the input image to segment the area around each maximum
    auto segmentation_areas = local_areas_voronoi_labels.clone();
    segmentation_areas << values::set_where_not(colors::grayscale32s::background(), local_max_areas_mask); // Do not consider areas outside of circle

    images::mask local_otsu_mask = semantic_convert<images::grayscale8u >(img) << binarize::otsu_per_component(segmentation_areas);

    // Slice the glomeruli by the voronoi borders
    local_areas_voronoi_labels << edge::laplacian();

    images::mask voronoi_borders(local_areas_voronoi_labels.get_image() != 0);
    local_otsu_mask << values::set_where(colors::mask::background(), voronoi_borders);

    // Blacklist phase: Get rid of everything we don't want
    // First calculate the 8-connected components (8 because this will help removing touching pixels)
    int local_otsu_mask_max_component_id = 0;
    images::grayscale32s local_otsu_mask_components = labeling::connected_components<colors::labels>(local_otsu_mask,
                                                                                                    local_otsu_mask_max_component_id,
                                                                                                    labeling::connectivity::connect_8);
    identity_recoloring_hashmap<colors::labels> blacklist;

    // Blacklist the components that border to the max circle and components that are too small
    segment_glomeruli_local_otsu_blacklist_by_contour(img,
                                                      blobs_maxima,
                                                      local_otsu_mask,
                                                      local_otsu_mask_components,
                                                      local_otsu_mask_max_component_id,
                                                      blacklist);

    local_otsu_mask_components << recoloring::recolor(blacklist);

//            images::mask local_otsu_mask_onlymaxima_noborder = toolbox2::to_mask(local_otsu_mask_components);
//            image_filter2::visualization::masks({local_otsu_mask,
//                                                 local_otsu_mask_onlymaxima_noborder,
//                                                 local_max_areas_mask,
//                                                 blobs_maxima,
//                                                 voronoi_borders}, img).show_and_wait("local otsu after non-maxima + no border blacklist vis");

    return toolbox::mask::from(local_otsu_mask_components);
}

void segmentation2d_local_otsu::create_parameters(misa_parameter_builder &t_parameters) {
    m_median_filter_size = t_parameters.create_algorithm_parameter<int>("median-filter-size", 3);
    m_glomeruli_min_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-min-rad", 15);
    m_glomeruli_max_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-max-rad", 65);
    m_cortex_segmentation_dilation_group_size = t_parameters.create_algorithm_parameter<double>("cortex-segmentation-dilation-group-size", 5);
    m_voronoi_cell_radius_border = t_parameters.create_algorithm_parameter<int>("voronoi-cell-radius-border", 5);
    m_isoperimetric_quotient_threshold = t_parameters.create_algorithm_parameter<double>("isoperimetric-quotient-threshold", 0.8);
}
