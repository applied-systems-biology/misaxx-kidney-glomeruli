//
// Created by rgerst on 11.09.18.
//


#pragma once

#include "segmentation2d_base.h"

namespace misaxx_kidney_glomeruli {
    struct segmentation2d_klingberg : public segmentation2d_base {

        using segmentation2d_base::segmentation2d_base;

        parameter<int> m_median_filter_size;
        parameter<double>  m_glomeruli_min_rad;
        parameter<double>  m_glomeruli_max_rad;
        parameter<double> m_threshold_percentile;
        parameter<double> m_threshold_factor;

        void work() override;

        void create_parameters(misaxx::misa_parameter_builder &t_parameters) override;
    };
}