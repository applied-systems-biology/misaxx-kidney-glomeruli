//
// Created by rgerst on 11.09.18.
//


#pragma once

#include "segmentation2d_base.h"

namespace misaxx_kidney_glomeruli {
    struct segmentation2d_klingberg : public segmentation2d_base {

        using segmentation2d_base::segmentation2d_base;

        int m_median_filter_size = from_algorithm_json_or<int>("median-filter-size", 3);
        double  m_glomeruli_min_rad = from_algorithm_json_or<double>("glomeruli-min-rad", 15);
        double  m_glomeruli_max_rad = from_algorithm_json_or<double>("glomeruli-max-rad", 65);
        double m_threshold_percentile = from_algorithm_json_or<double>("threshold-percentile", 75);
        double m_threshold_factor = from_algorithm_json_or<double>("threshold-factor", 1.5);

        void misa_work() override;
    };
}