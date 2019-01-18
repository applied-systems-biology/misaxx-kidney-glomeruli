//
// Created by rgerst on 13.09.18.
//


#pragma once

#include "segmentation2d_base.h"
#include <misaxx-imaging/coixx/toolbox/toolbox_blob.h>
#include <misaxx-imaging/coixx/toolbox/toolbox_localminmax.h>
#include <misaxx-imaging/coixx/toolbox/toolbox_binarize_componentotsu.h>
#include <misaxx-imaging/coixx/recoloring_map.h>

namespace misaxx_kidney_glomeruli {
    struct segmentation2d_local_otsu : public segmentation2d_base {

        parameter<int> m_median_filter_size;
        parameter<double>  m_glomeruli_min_rad;
        parameter<double>  m_glomeruli_max_rad;
        parameter<double> m_cortex_segmentation_dilation_group_size;
        parameter<int> m_voronoi_cell_radius_border;
        parameter<double> m_isoperimetric_quotient_threshold;

        using segmentation2d_base::segmentation2d_base;

        void work() override;

        void create_parameters(misaxx::misa_parameter_builder &t_parameters) override;

    protected:

        /**
         * Extracts blobs using LoG
         * @param img
         * @param borders
         * @return
         */
        coixx::images::grayscale_float extract_blobs_log(const coixx::images::grayscale_float &img);

        /**
         * Finds the local maxima
         * @param blobs
         * @param img
         * @param tissue_mask
         * @return
         */
        coixx::images::mask extract_maxima(coixx::images::grayscale_float blobs,
                                    const coixx::images::mask &tissue_mask,
                                    const double t_radius);

        /**
        * Given a mask of local maxima, only select the maxima that are within a larger radius the input
        * @param t_maxima
        * @param t_radius
        * @return
        */
        void restrict_maxima(coixx::images::mask &t_maxima, const coixx::images::grayscale_float &img, const double t_radius);

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
        coixx::images::mask find_cortex_otsu_distance_and_dilation(const coixx::images::mask &tissue_mask, const coixx::images::grayscale_float &blobs);

        /**
         * Uses OpenCV's contours feature to remove non-glomeruli
         * @param blobs_maxima_mask
         * @param local_otsu_mask_components
         * @param blacklist
         */
        void segment_glomeruli_local_otsu_blacklist_by_contour(const coixx::images::grayscale_float &img,
                                                               const coixx::images::mask &blobs_maxima_mask,
                                                               const coixx::images::mask &local_otsu_mask,
                                                               const coixx::images::grayscale32s &local_otsu_mask_components,
                                                               const int local_otsu_mask_max_component_id,
                                                               coixx::mutable_recoloring_map<coixx::colors::labels> &blacklist);

        /**
         * Uses a local Otsu approach to segment the glomeruli
         * @return
         */
        coixx::images::mask segment_glomeruli_local_otsu(const coixx::images::mask &blobs_maxima, const coixx::images::mask &cortex_mask, const coixx::images::grayscale_float &img);
    };
}