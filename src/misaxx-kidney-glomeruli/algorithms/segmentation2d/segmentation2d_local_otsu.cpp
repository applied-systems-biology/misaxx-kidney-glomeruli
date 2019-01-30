//
// Created by rgerst on 17.12.18.
//

#include "segmentation2d_local_otsu.h"
#include <cv-toolbox/ReadableBMatTypes.h>
#include <cv-toolbox/toolbox/toolbox_blob.h>
#include <cv-toolbox/toolbox/toolbox_statistics.h>
#include <cv-toolbox/toolbox/toolbox_values.h>
#include <cv-toolbox/toolbox/toolbox_semantic_convert.h>
#include <cv-toolbox/toolbox/toolbox_binarize.h>
#include <cv-toolbox/toolbox/toolbox_morph.h>
#include <cv-toolbox/toolbox/toolbox_resize.h>
#include <cv-toolbox/toolbox/toolbox_edge.h>
#include <cv-toolbox/structuring_element.h>
#include <cv-toolbox/recoloring_map.h>
#include <cv-toolbox/toolbox/toolbox_recoloring.h>
#include <cv-toolbox/toolbox/toolbox_normalize.h>

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

namespace {
    cv::images::grayscale32f extract_blobs_log(const cv::images::grayscale32f &img,
            double glomeruli_min_rad, double glomeruli_max_rad, double voxel_xy) {

        const double glomeruli_min_rad_sigma = (glomeruli_min_rad / voxel_xy) / sqrt(2);
        const double glomeruli_max_rad_sigma = (glomeruli_max_rad / voxel_xy) / sqrt(2);
        const double avg_rad_sigma = (glomeruli_min_rad_sigma + glomeruli_max_rad_sigma) / 2;

        cv::images::grayscale32f log_response = img.clone();
        cv::toolbox::laplacian_of_gaussian(log_response, avg_rad_sigma);

        // Manual postprocessing for increased speed:
        // Takes only the negative response and normalizes against the lowest response (where the blobs are)
        float min_response = static_cast<float>(cv::toolbox::statistics::min_max_loc(log_response).min.value);
        for(int y = 0; y < img.rows; ++y) {
            float *row = log_response[y];
            for(int x = 0; x < img.cols; ++x) {
                row[x] = std::min(0.0f, row[x]) / min_response;
            }
        }

        return log_response;
    }

    cv::images::mask extract_maxima(const cv::images::grayscale32f &blobs, const cv::images::mask &tissue_mask,
                                                                  const double t_radius) {

        // Find the local maxima
        cv::images::mask blobs_mask = cv::toolbox::exclusive_local_maxima(blobs, static_cast<int>(2 * t_radius));
        cv::toolbox::set_where_not<uchar>(blobs_mask, tissue_mask, 0);

        return blobs_mask;
    }

    void restrict_maxima(cv::images::mask &t_maxima, const cv::images::grayscale32f &img,
                                               const double t_radius) {

        // Get all maxima
        std::vector<cv::pixel<float>> pixels;

        for(int y = 0; y < t_maxima.rows; ++y) {
            const uchar *row = t_maxima[y];
            const float *row_values = img[y];
            for(int x = 0; x < t_maxima.cols; ++x ) {
                if(row[x] > 0) {
                    pixels.emplace_back(cv::pixel<float>(x, y, row_values[x]));
                }
            }
        }

        const double r_sq = pow(2 * t_radius, 2);

        // For each pixel, check if we have another pixel with same or larger value within radius
        // If this is the case, we drop the pixel
        for(size_t i = 0; i < pixels.size(); ++i) {
            const auto px = pixels[i];
            for(size_t j = 0; j < pixels.size(); ++j) {
                if(i != j) {
                    const auto px2 = pixels[j];
                    const double l = pow(px.location.x - px2.location.x, 2) + pow(px.location.y - px2.location.y, 2);
                    if(l <= r_sq && px2.value >= px.value) {
                        t_maxima.at<uchar>(px.location) = 0; // Get rid of the maximum
                        break;
                    }
                }
            }
        }
    }

    cv::images::mask find_cortex_otsu_distance_and_dilation(const cv::images::mask &tissue_mask,
                                                                      const cv::images::grayscale32f &blobs,
                                                                      double voxel_xy,
                                                                      double glomeruli_max_rad,
                                                                      double cortex_segmentation_dilation_group_size) {

        const int glomeruli_max_diameter = static_cast<int>((glomeruli_max_rad / voxel_xy) * 2);

        cv::images::grayscale32f cortex_dist = cv::images::grayscale32f::allocate(blobs.size());
        cv::distanceTransform(tissue_mask, cortex_dist, cv::DIST_L2, 3);

        // Apply otsu
        cv::images::grayscale8u blobs_thresholded = cv::toolbox::semantic_convert::to_grayscale8u(blobs);
        cv::toolbox::otsu_where(blobs_thresholded, tissue_mask);

        // Dilation based method
        int dilate_diameter = static_cast<int>(glomeruli_max_diameter * 0.1 * cortex_segmentation_dilation_group_size);
        cv::images::grayscale8u blobs_thresholded_small =
                cv::toolbox::resize(blobs_thresholded, 0.1, cv::toolbox::resize_interpolation::nearest);
        cv::toolbox::morph::dilate(blobs_thresholded_small, cv::structuring_element::ellipse(dilate_diameter));

        cv::images::mask cortex_by_dilation = cv::toolbox::resize(blobs_thresholded_small,
                                                 blobs_thresholded.size(),
                                                 cv::toolbox::resize_interpolation::nearest);
        cv::toolbox::set_where_not<uchar>(cortex_by_dilation, tissue_mask, 0);
        cv::toolbox::otsu(blobs_thresholded);

        // Narrow down to glomeruli regions & tissue
        auto cortex_dist_candidates = cortex_dist.clone();
        cv::toolbox::set_where_not<float>(cortex_dist_candidates, tissue_mask, 0.0f);
        cv::toolbox::set_where_not<float>(cortex_dist_candidates, blobs_thresholded, 0.0f);

        // Options: Maximum (works well in test images, BUT outliers disturb it!)
        // Alternative: 2 * mean (seems to work)
//            float dist_max = toolbox(cortex_dist_candidates).max();
        cv::images::mask cortex_dist_candidates_nonzero { cortex_dist_candidates > 0 };
        double dist = cv::mean(cortex_dist_candidates, cortex_dist_candidates_nonzero)[0] * 2;

        cv::images::mask cortex { cortex_dist <= dist };

        // Combine both
        cv::toolbox::set_where<uchar>(cortex, cortex_by_dilation, 255);
        cv::toolbox::set_where_not<uchar>(cortex, tissue_mask, 0);

        return cortex;
    }

    void segment_glomeruli_local_otsu_blacklist_by_contour(const cv::images::grayscale32f &,
                                                             const cv::images::mask &blobs_maxima_mask,
                                                             const cv::images::mask &local_otsu_mask,
                                                             const cv::images::grayscale32s &local_otsu_mask_components,
                                                             const int local_otsu_mask_max_component_id,
                                                             cv::mutable_recoloring_map<int> &blacklist,
                                                             const double voxel_xy,
                                                             const double glomeruli_min_rad_,
                                                             const double glomeruli_max_rad_,
                                                             const double isoperimetric_quotient_threshold) {

        // Detect the contours
        std::vector<std::vector<cv::Point>> contours;
        cv::findContours(local_otsu_mask, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

        // Find the center for each label
        std::vector<cv::Point> label_centers;
        label_centers.resize(static_cast<size_t>(local_otsu_mask_max_component_id) + 1, cv::Point(-1, -1));

        for(int y = 0; y < local_otsu_mask_components.rows; ++y) {

            const uchar* row_maxima = blobs_maxima_mask[y];
            const int* row_components = local_otsu_mask_components[y];

            for(int x = 0; x < local_otsu_mask_components.cols; ++x) {
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
        const int glomeruli_min_rad = static_cast<int>(glomeruli_min_rad_ / voxel_xy);
        const int glomeruli_max_rad = static_cast<int>(glomeruli_max_rad_ / voxel_xy);

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
                    blacklist.set_recolor(label, 0);
//                        std::cout << "Blacklist " << label << " failed radius requirements" << std::endl;
                    continue;
                }

                // Eccentricity check
                if(max_rad / min_rad >= 2) {
                    blacklist.set_recolor(label, 0);
//                        std::cout << "Blacklist " << label << " failed eccentricity requirements" << std::endl;
                    continue;
                }

                // Countour roundness check
                const double isoperimetric_quotient = (4 * M_PI * cv::contourArea(contour)) / pow(cv::arcLength(contour, true), 2);
                const double fitted_ellipse_area = M_PI * fitted_ellipse_r1 * fitted_ellipse_r2;
                const double fitted_ellipse_perimeter = 2 * M_PI * sqrt((pow(fitted_ellipse_r1, 2) + pow(fitted_ellipse_r2, 2)) / 2.0);
                const double fitted_isoperimetric_quotient = (4 * M_PI * fitted_ellipse_area) / pow(fitted_ellipse_perimeter, 2);

//                    std::cout << "Q = " << isoperimetric_quotient << ", expected Q = " << fitted_isoperimetric_quotient << std::endl;

                if(isoperimetric_quotient < isoperimetric_quotient_threshold * fitted_isoperimetric_quotient || isoperimetric_quotient > 1) {
                    blacklist.set_recolor(label, 0);
//                        std::cout << "Blacklist " << label << " failed isoperimetric quotient requirements" << std::endl;
                    continue;
                }

            }
            else {
//                    std::cout << "Blacklist " << label << " has no contour" << std::endl;
                blacklist.set_recolor(label, 0);
            }
        }

    }

    cv::images::mask segment_glomeruli_local_otsu(const cv::images::mask &blobs_maxima,
            const cv::images::mask &cortex_mask,
            const cv::images::grayscale32f &img,
            const double voxel_xy,
            const double glomeruli_max_rad_,
            const double glomeruli_min_rad_,
            const int voronoi_cell_radius_border,
            const double isoperimetric_quotient_threshold) {

        const int glomeruli_max_diameter = static_cast<int>((glomeruli_max_rad_ / voxel_xy) * 2);
        const int glomeruli_search_diameter = glomeruli_max_diameter + voronoi_cell_radius_border;

        // Create an area around each maximum
        cv::images::mask local_max_areas_mask = blobs_maxima.clone();
        cv::toolbox::morph::dilate(local_max_areas_mask, cv::structuring_element::ellipse(glomeruli_search_diameter));

        cv::toolbox::set_where_not<uchar>(local_max_areas_mask, cortex_mask, 0); // Remove regions outside cortex

        // Generate voronoi partitioning
        cv::images::grayscale32f local_areas_voronoi_distance(img.size(), 0);
        cv::images::grayscale8u local_areas_voronoi_centers = blobs_maxima.clone();
        cv::toolbox::invert(local_areas_voronoi_centers);
        cv::images::grayscale32s local_areas_voronoi_labels(img.size(), 0);
        cv::distanceTransform(local_areas_voronoi_centers,
                              local_areas_voronoi_distance,
                              local_areas_voronoi_labels,
                              cv::DIST_L2,
                              3);

        // Delete partitioning outside cortex
        cv::toolbox::set_where_not<int>(local_areas_voronoi_labels, cortex_mask, 0);

        // Apply per-component localized Otsu on the input image to segment the area around each maximum
        cv::images::labels segmentation_areas = local_areas_voronoi_labels.clone();
        cv::toolbox::set_where_not(segmentation_areas, local_max_areas_mask, 0); // Do not consider areas outside of circle

        cv::images::mask local_otsu_mask = cv::toolbox::semantic_convert::to_grayscale8u(img);
        cv::toolbox::otsu_per_component(local_otsu_mask, segmentation_areas);

        // Slice the glomeruli by the voronoi borders
        cv::toolbox::laplacian<int>(local_areas_voronoi_labels);

        cv::images::mask voronoi_borders(local_areas_voronoi_labels != 0);
        cv::toolbox::set_where<uchar>(local_otsu_mask, voronoi_borders, 0);

        // Blacklist phase: Get rid of everything we don't want
        // First calculate the 8-connected components (8 because this will help removing touching pixels)
        int local_otsu_mask_max_component_id = 0;
        cv::images::grayscale32s local_otsu_mask_components;
        cv::connectedComponents(local_otsu_mask, local_otsu_mask_components, 8, CV_32S);

        cv::identity_recoloring_hashmap<int> blacklist;

        // Blacklist the components that border to the max circle and components that are too small
        segment_glomeruli_local_otsu_blacklist_by_contour(img,
                                                          blobs_maxima,
                                                          local_otsu_mask,
                                                          local_otsu_mask_components,
                                                          local_otsu_mask_max_component_id,
                                                          blacklist, voxel_xy,
                                                          glomeruli_min_rad_,
                                                          glomeruli_max_rad_,
                                                          isoperimetric_quotient_threshold);

        cv::toolbox::recolor(local_otsu_mask_components, blacklist);


//            images::mask local_otsu_mask_onlymaxima_noborder = toolbox2::to_mask(local_otsu_mask_components);
//            image_filter2::visualization::masks({local_otsu_mask,
//                                                 local_otsu_mask_onlymaxima_noborder,
//                                                 local_max_areas_mask,
//                                                 blobs_maxima,
//                                                 voronoi_borders}, img).show_and_wait("local otsu after non-maxima + no border blacklist vis");

        return cv::toolbox::to_mask(local_otsu_mask_components);
    }
}

void segmentation2d_local_otsu::work() {

    auto module = get_module_as<module_interface>();

    cv::images::mask tissue_mask = cv::toolbox::semantic_convert::to_grayscale8u(m_input_tissue.clone());

    if(cv::countNonZero(tissue_mask) == 0) { //INFO: Not inverted yet
        // Instead save a black image
        cv::images::mask img { tissue_mask.size(), 0 };
        m_output_segmented2d.write(std::move(img));
        return;
    }

    const double voxel_xy = module->m_voxel_size.get_size_xy().get_value();

    cv::images::grayscale32f img = cv::toolbox::semantic_convert::to_grayscale32f(m_input_autofluoresence.clone());
    cv::images::grayscale32f img_original = img.clone();

    // Smooth & normalize
    cv::medianBlur(img, img.buffer(), m_median_filter_size.query()); img.swap();
    cv::toolbox::normalize::by_max(img);

    // Extracts the blobs
    cv::images::grayscale32f blobs = extract_blobs_log(img,
            m_glomeruli_min_rad.query(),
            m_glomeruli_max_rad.query(),
            voxel_xy);
    cv::images::mask cortex_mask = find_cortex_otsu_distance_and_dilation(tissue_mask,
            blobs,
            voxel_xy,
            m_glomeruli_max_rad.query(),
            m_cortex_segmentation_dilation_group_size.query());

    cv::toolbox::set_where_not<float>(img, cortex_mask, 0);

    // Run multiple iterations of local Otsu and merge them if there are no conflicts
    cv::images::mask merged_mask { img.size(), 0 };

    // Local maxima selection is expensive (Due to dilation).
    // Instead use a two-step approach that only requires the dilation with the small selection
    cv::images::mask blobs_all_maxima = extract_maxima(blobs.clone(), tissue_mask, m_glomeruli_min_rad.query() / voxel_xy);
    cv::toolbox::set_where_not<uchar>(blobs_all_maxima, cortex_mask, 0);

    // We first try to segment large glomeruli to prevent oversegmentation (in combination with deleting already segmented maxima)
    for(const double radius_microns : { m_glomeruli_max_rad.query(), m_glomeruli_min_rad.query() }) {
        const double radius = radius_microns / voxel_xy;

        cv::images::mask blobs_maxima = blobs_all_maxima.clone();
        if(radius_microns != m_glomeruli_min_rad.query()) {
            restrict_maxima(blobs_maxima, img, radius);
        }

        cv::images::grayscale32s blobs_maxima_components;
        cv::connectedComponents(blobs_maxima, blobs_maxima_components, 8, CV_32S);

        cv::images::mask final_mask = segment_glomeruli_local_otsu(blobs_maxima,
                cortex_mask,
                img,
                voxel_xy,
                m_glomeruli_max_rad.query(),
                m_glomeruli_min_rad.query(),
                m_voronoi_cell_radius_border.query(),
                m_isoperimetric_quotient_threshold.query());

        // Delete only good positions from the list of maxima to be analyzed.
        // The reason behind this is that a larger search radius can cover two adjacent glomeruli where
        // only one maximum is seen as relevant. This will lead to bad detection. But if the maximum is already
        // deleted, the small scale algorithm won't be able to detect the right glomeruli.
        cv::toolbox::set_where<uchar>(blobs_all_maxima, final_mask, 0);

        cv::bitwise_or(merged_mask, final_mask, merged_mask.buffer());
        merged_mask.swap();
    }

//            image_filter2::visualization::mask(merged_mask, img).show_and_wait("final mask");
    m_output_segmented2d.write(std::move(merged_mask));
}

void segmentation2d_local_otsu::create_parameters(misa_parameter_builder &t_parameters) {
    segmentation2d_base::create_parameters(t_parameters);
    m_median_filter_size = t_parameters.create_algorithm_parameter<int>("median-filter-size", 3);
    m_glomeruli_min_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-min-rad", 15);
    m_glomeruli_max_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-max-rad", 65);
    m_cortex_segmentation_dilation_group_size = t_parameters.create_algorithm_parameter<double>("cortex-segmentation-dilation-group-size", 5);
    m_voronoi_cell_radius_border = t_parameters.create_algorithm_parameter<int>("voronoi-cell-radius-border", 5);
    m_isoperimetric_quotient_threshold = t_parameters.create_algorithm_parameter<double>("isoperimetric-quotient-threshold", 0.8);
}
