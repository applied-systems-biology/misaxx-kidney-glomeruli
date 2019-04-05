/**
 * Copyright by Ruman Gerst
 * Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
 * https://www.leibniz-hki.de/en/applied-systems-biology.html
 * HKI-Center for Systems Biology of Infection
 * Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Insitute (HKI)
 * Adolf-Reichwein-Straße 23, 07745 Jena, Germany
 *
 * This code is licensed under BSD 2-Clause
 * See the LICENSE file provided with this code for the full license.
 */

#include "segmentation2d_klingberg.h"

using namespace misaxx;
using namespace misaxx_kidney_glomeruli;

namespace cv::images {
    using grayscale8u = cv::Mat1b;
    using grayscale32f = cv::Mat1f;
    using mask = cv::Mat1b;
}

namespace {

    template<typename T>
    std::vector<T> get_sorted_pixels(const cv::Mat_<T> &img) {
        std::vector<T> pixels;
        pixels.reserve(img.rows * img.cols);
        for(int y = 0; y < img.rows; ++y) {
            const auto *row = img[y];
            for(int x = 0; x < img.cols; ++x) {
                pixels.push_back(row[x]);
            }
        }
        std::sort(pixels.begin(), pixels.end());
        return pixels;
    }

    template<typename T>
    std::vector<T> get_percentiles(const std::vector<T> &pixels, const std::vector<double> &percentiles) {
        std::vector<T> result;
        for(double percentile : percentiles) {
            size_t rank = static_cast<size_t>(std::ceil(percentile / 100.0 * pixels.size()));
            result.push_back(pixels.at(rank));
        }
        return result;
    }

    double get_max_value(const cv::Mat &img) {
        double max;
        cv::minMaxLoc(img, nullptr, &max);
        return max;
    }

    void normalize_by_max(cv::images::mask &img) {
        const double max = get_max_value(img);
        for(int y = 0; y < img.rows; ++y) {
            auto *row = img[y];
            for(int x = 0; x < img.cols; ++x) {
                row[x] = static_cast<uchar>(row[x] * 255.0 / max);
            }
        }
    }

    void normalize_by_max(cv::images::grayscale32f &img) {
        const double max = get_max_value(img);
        for(int y = 0; y < img.rows; ++y) {
            auto *row = img[y];
            for(int x = 0; x < img.cols; ++x) {
                row[x] = static_cast<float>(row[x] / max);
            }
        }
    }

    cv::images::grayscale32f get_as_grayscale_float_copy(const cv::Mat &img) {
        if(img.type() == CV_32F) {
            return img.clone();
        }
        else if(img.type() == CV_64F) {
            cv::images::grayscale32f result;
            img.convertTo(result, CV_32F, 1);
            return result;
        }
        else if(img.type() == CV_8U) {
            cv::images::grayscale32f result;
            img.convertTo(result, CV_32F, 1.0 / 255.0);
            return result;
        }
        else if(img.type() == CV_16U) {
            cv::images::grayscale32f result;
            img.convertTo(result, CV_32F, 1.0 / std::numeric_limits<ushort>::max());
            return result;
        }
        else {
            throw std::runtime_error("Unsupported image depth: " + std::to_string(img.type()));
        }
    }

    cv::images::grayscale8u get_preprocessed_image(const misaxx::ome::misa_ome_plane &plane, int median_filter_size) {
        cv::images::grayscale32f img = get_as_grayscale_float_copy(plane.access_readonly().get());

        // Initial median filtering + normalization
        cv::medianBlur(img.clone(), img, median_filter_size);
        normalize_by_max(img);

        cv::images::grayscale8u img8u;
        img.convertTo(img8u, CV_8U, 255.0);
        return img8u;
    }

}

void segmentation2d_klingberg::work() {

    auto pam = is_parallelizeable_parameter;

    auto module = get_module_as<module_interface>();

    cv::images::mask img_non_tissue_mask { m_input_tissue.clone() };

    if(cv::countNonZero(img_non_tissue_mask) == 0) { //INFO: Not inverted yet
        // Instead save a black image
        m_output_segmented2d.write(cv::images::mask(img_non_tissue_mask.size(), 0));
        return;
    }

    img_non_tissue_mask = 255 - img_non_tissue_mask; // Now we mark all non-tissue

    // Get the preprocessed image
    cv::images::grayscale8u img8u = get_preprocessed_image(m_input_autofluoresence, m_median_filter_size.query());

    // Generated parameters
    const double voxel_xy = module->m_voxel_size.get_size_xy().get_value();
    int glomeruli_max_morph_disk_radius = static_cast<int>(m_glomeruli_max_rad.query() / voxel_xy);
    int glomeruli_min_morph_disk_radius = static_cast<int>((m_glomeruli_min_rad.query() / 2.0) / voxel_xy);

    // Morphological operation (opening)
    // Corresponds to only allowing objects > disk_size to be included
    // Also subtract the morph result from the initial to remove uneven background + normalize
    cv::morphologyEx(img8u.clone(), img8u, cv::MORPH_TOPHAT, cv::getStructuringElement(cv::MORPH_ELLIPSE,
            cv::Size(glomeruli_max_morph_disk_radius * 2 + 1, glomeruli_max_morph_disk_radius * 2 + 1)));
    normalize_by_max(img8u);

    // We are first extracting tissue data
    cv::images::grayscale8u img_tissue = img8u.clone();

    // Get rid of non-tissue
    img_tissue.setTo(0, img_non_tissue_mask);

    // Only select glomeruli if the threshold is higher than 75-percentile of kidney tissue
    double percentile_tissue = get_percentiles(get_sorted_pixels(img8u), { m_threshold_percentile.query() })[0];

    //////////////
    // Now working in uint8
    //////////////

    // Threshold the main image
    double otsu_threshold = cv::threshold(img8u.clone(), img8u, 0, 255, cv::THRESH_OTSU);

//    std::cout << "Otsu: " << std::to_string(otsu_threshold) << " Percentile: " << std::to_string(percentile_tissue) << std::endl;

    if(otsu_threshold > percentile_tissue * m_threshold_factor.query() ) {

        // Get rid of non-tissue
        img8u.setTo(0, img_non_tissue_mask);

        // Morphological operation (object should have min. radius)
        cv::morphologyEx(img8u.clone(), img8u, cv::MORPH_OPEN, cv::getStructuringElement(cv::MORPH_ELLIPSE,
                cv::Size(glomeruli_min_morph_disk_radius * 2 + 1, glomeruli_min_morph_disk_radius * 2 + 1)));
    }
    else {
        img8u = 0;
    }

    // Save the mask
    m_output_segmented2d.write(std::move(img8u));
}

void segmentation2d_klingberg::create_parameters(misa_parameter_builder &t_parameters) {
    segmentation2d_base::create_parameters(t_parameters);
    m_median_filter_size = t_parameters.create_algorithm_parameter<int>("median-filter-size", 3);
    m_glomeruli_min_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-min-rad", 15);
    m_glomeruli_max_rad = t_parameters.create_algorithm_parameter<double>("glomeruli-max-rad", 65);
    m_threshold_percentile = t_parameters.create_algorithm_parameter<double>("threshold-percentile", 75);
    m_threshold_factor = t_parameters.create_algorithm_parameter<double>("threshold-factor", 1.5);
}
