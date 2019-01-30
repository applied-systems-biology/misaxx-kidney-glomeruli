//
// Created by rgerst on 13.09.18.
//


#pragma once

#include "segmentation2d_base.h"

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

    };
}