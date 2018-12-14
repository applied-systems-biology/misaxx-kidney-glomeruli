//
// Created by rgerst on 13.09.18.
//


#pragma once

#include "segmentation2d_base.h"
#include <coixx/toolbox/toolbox_blob.h>
#include <coixx/toolbox/toolbox_localminmax.h>
#include <coixx/toolbox/toolbox_binarize_componentotsu.h>
#include <misaxx/vdata/object3d_voxel_size.h>

namespace misaxx_kidney_glomeruli {
    struct segmentation2d_local_otsu : public segmentation2d_base {

        int m_median_filter_size = from_algorithm_json_or<int>("median-filter-size", 3);
        double  m_glomeruli_min_rad = from_algorithm_json_or<double>("glomeruli-min-rad", 15);
        double  m_glomeruli_max_rad = from_algorithm_json_or<double>("glomeruli-max-rad", 65);
        double m_cortex_segmentation_dilation_group_size = from_algorithm_json_or<double>("cortex-segmentation-dilation-group-size", 5);
        int m_voronoi_cell_radius_border = from_algorithm_json_or<int>("voronoi-cell-radius-border", 5);
        double m_isoperimetric_quotient_threshold = from_algorithm_json_or<double>("isoperimetric-quotient-threshold", 0.8);
        object3d_voxel_size m_voxel_size = from_parameter<object3d_voxel_size>();

        using segmentation2d_base_base;

        void misa_work() override {
           using namespace coixx::toolbox;

            images::mask tissue_mask = m_input_tissue->load();

            if(statistics::is_black(tissue_mask)) { //INFO: Not inverted yet
                // Instead save a black image
                images::mask img(tissue_mask.get_size(), colors::mask::background());
                m_output_segmented2d->save(img);
                return;
            }

            images::grayscale_float img = m_input_autofluoresence->auto_load();
            images::grayscale_float img_original;

            // Smooth & normalize
            img << values::backup(img_original) << blur::median(m_median_filter_size) << normalize::by_max();

            // Extracts the blobs
            images::grayscale_float blobs = extract_blobs_log(img);

            images::mask cortex_mask = find_cortex_otsu_distance_and_dilation(tissue_mask, blobs);
//            image_filter2::visualization::mask(cortex_mask, img).show_and_wait("cortex visualization");


            img << values::set_where_not(colors::grayscale_float::black(), cortex_mask);

            // Run multiple iterations of local Otsu and merge them if there are no conflicts
            images::mask merged_mask(img.get_size(), colors::mask::background());

            // Local maxima selection is expensive (Due to dilation).
            // Instead use a two-step approach that only requires the dilation with the small selection
            images::mask blobs_all_maxima = extract_maxima(blobs.clone(), img, tissue_mask, m_glomeruli_min_rad / m_voxel_size.xz());
            blobs_all_maxima << values::set_where_not(colors::mask::background(), cortex_mask);

            // We first try to segment large glomeruli to prevent oversegmentation (in combination with deleting already segmented maxima)
            for(const double radius_microns : { m_glomeruli_max_rad, m_glomeruli_min_rad }) {
                const double radius = radius_microns / m_voxel_size.xz();

                images::mask blobs_maxima = blobs_all_maxima.clone();
                if(radius_microns != m_glomeruli_min_rad) {
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
            m_output_segmented2d->save(merged_mask);
        }

    protected:

        /**
         * Extracts blobs using LoG
         * @param img
         * @param borders
         * @return
         */
        images::grayscale_float extract_blobs_log(const images::grayscale_float &img) {

            using namespace coixx::toolbox;

            const double glomeruli_min_rad_sigma = (m_glomeruli_min_rad / m_voxel_size.xz()) / sqrt(2);
            const double glomeruli_max_rad_sigma = (m_glomeruli_max_rad / m_voxel_size.xz()) / sqrt(2);
            const double avg_rad_sigma = (glomeruli_min_rad_sigma + glomeruli_max_rad_sigma) / 2;

            images::grayscale_float log_response = img.clone() << blob::laplacian_of_gaussian(avg_rad_sigma);

            // Manual postprocessing for increased speed:
            // Takes only the negative response and normalizes against the lowest response (where the blobs are)
            float min_response = statistics::min(log_response);
            for(int y = 0; y < img.get_image().rows; ++y) {
                float *row = log_response.raw_row_ptr(y);
                for(int x = 0; x < img.get_image().cols; ++x) {
                    row[x] = std::min(0.0f, row[x]) / min_response;
                }
            }

            return log_response;
        }

        /**
         * Finds the local maxima
         * @param blobs
         * @param img
         * @param tissue_mask
         * @return
         */
        images::mask extract_maxima(images::grayscale_float blobs,
                                    const images::grayscale_float &img,
                                    const images::mask &tissue_mask,
                                    const double t_radius) {

            using namespace coixx::toolbox;

            // Find the local maxima
            images::mask blobs_mask = localminmax::local_exclusive_max_morph(blobs, static_cast<int>(2 * t_radius));
            blobs_mask << values::set_where_not(colors::mask::black(), tissue_mask);

            return blobs_mask;
        }

        /**
        * Given a mask of local maxima, only select the maxima that are within a larger radius the input
        * @param t_maxima
        * @param t_radius
        * @return
        */
        void restrict_maxima(images::mask &t_maxima, const images::grayscale_float &img, const double t_radius) {

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

        /**
         * Finds cortex by combining two methods:
         * Using Otsu to find "true" glomeruli and measuring the mean distance to the tissue border.
         * This distance is used to calculate the border that contains the glomeruli.
         * The second method dilates the "true" glomeruli positions morphologically and creates an area of
         * glomeruli.
         *
         * This increases the chance of finding strong glomeruli that are more inside of the kidney.
         *
         * @param tissue_mask
         * @param blobs
         * @return
         */
        images::mask find_cortex_otsu_distance_and_dilation(const images::mask &tissue_mask, const images::grayscale_float &blobs) {

            using namespace coixx::toolbox;

            const int glomeruli_max_diameter = static_cast<int>((m_glomeruli_max_rad / m_voxel_size.xz()) * 2);

            images::grayscale_float cortex_dist(blobs.get_size(), colors::grayscale_float::black());
            cv::distanceTransform(tissue_mask.get_image(), cortex_dist.get_image(), cv::DIST_L2, 3);

            // Apply otsu
            images::grayscale8u blobs_thresholded = semantic_convert<images::grayscale8u >(blobs) << binarize::otsu_where(tissue_mask);

            // Dilation based method
            int dilate_diameter = static_cast<int>(glomeruli_max_diameter * 0.1 * m_cortex_segmentation_dilation_group_size);
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

        /**
         * Uses OpenCV's contours feature to remove non-glomeruli
         * @param blobs_maxima_mask
         * @param local_otsu_mask_components
         * @param blacklist
         */
        void segment_glomeruli_local_otsu_blacklist_by_contour(const images::grayscale_float &img,
                                                               const images::mask &blobs_maxima_mask,
                                                               const images::mask &local_otsu_mask,
                                                               const images::grayscale32s &local_otsu_mask_components,
                                                               const int local_otsu_mask_max_component_id,
                                                               coixx::toolbox::recoloring::recoloring_map<images::labels> &blacklist) {

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

                const uchar * row_maxima = blobs_maxima_mask.raw_row_ptr(y);
                const int * row_components = local_otsu_mask_components.raw_row_ptr(y);

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
            const int glomeruli_min_rad = static_cast<int>(m_glomeruli_min_rad / m_voxel_size.xz());
            const int glomeruli_max_rad = static_cast<int>(m_glomeruli_max_rad / m_voxel_size.xz());

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
                        blacklist.define(static_cast<int>(label)) = 0;
//                        std::cout << "Blacklist " << label << " failed radius requirements" << std::endl;
                        continue;
                    }

                    // Eccentricity check
                    if(max_rad / min_rad >= 2) {
                        blacklist.define(static_cast<int>(label)) = 0;
//                        std::cout << "Blacklist " << label << " failed eccentricity requirements" << std::endl;
                        continue;
                    }

                    // Countour roundness check
                    const double isoperimetric_quotient = (4 * M_PI * cv::contourArea(contour)) / pow(cv::arcLength(contour, true), 2);
                    const double fitted_ellipse_area = M_PI * fitted_ellipse_r1 * fitted_ellipse_r2;
                    const double fitted_ellipse_perimeter = 2 * M_PI * sqrt((pow(fitted_ellipse_r1, 2) + pow(fitted_ellipse_r2, 2)) / 2.0);
                    const double fitted_isoperimetric_quotient = (4 * M_PI * fitted_ellipse_area) / pow(fitted_ellipse_perimeter, 2);

//                    std::cout << "Q = " << isoperimetric_quotient << ", expected Q = " << fitted_isoperimetric_quotient << std::endl;

                    if(isoperimetric_quotient < m_isoperimetric_quotient_threshold * fitted_isoperimetric_quotient || isoperimetric_quotient > 1) {
                        blacklist.define(static_cast<int>(label)) = 0;
//                        std::cout << "Blacklist " << label << " failed isoperimetric quotient requirements" << std::endl;
                        continue;
                    }

                }
                else {
//                    std::cout << "Blacklist " << label << " has no contour" << std::endl;
                    blacklist.define(static_cast<int>(label)) = 0;
                }
            }

        }

        /**
         * Uses a local Otsu approach to segment the glomeruli
         * @return
         */
        images::mask segment_glomeruli_local_otsu(const images::mask &blobs_maxima, const images::mask &cortex_mask, const images::grayscale_float &img) {

            using namespace coixx::toolbox;

            const int glomeruli_max_diameter = static_cast<int>((m_glomeruli_max_rad / m_voxel_size.xz()) * 2);
            const int glomeruli_search_diameter = glomeruli_max_diameter + m_voronoi_cell_radius_border;

            // Create an area around each maximum
            images::mask local_max_areas_mask = blobs_maxima.clone() << morph::dilate(structuring_element::ellipse(glomeruli_search_diameter));
            local_max_areas_mask << values::set_where_not(colors::mask::background(), cortex_mask); // Remove regions outside cortex

            // Generate voronoi partitioning
            images::grayscale_float local_areas_voronoi_distance(img.get_size(), colors::grayscale_float::black());
            images::grayscale8u local_areas_voronoi_centers = blobs_maxima.clone() << values::invert();
            images::grayscale32s local_areas_voronoi_labels(img.get_size());
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
            images::grayscale32s local_otsu_mask_components = labeling::connected_components<color_grayscale_int32>(local_otsu_mask,
                                                                                                                    local_otsu_mask_max_component_id,
                                                                                                                    labeling::connectivity::connect_8);
            recoloring::recoloring_map<images::labels> blacklist;
            blacklist.data().reserve(local_otsu_mask_max_component_id);

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

            return  mask::from(local_otsu_mask_components);
        }
    };
}