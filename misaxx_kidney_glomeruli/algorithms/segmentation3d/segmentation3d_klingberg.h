//
// Created by rgerst on 12.09.18.
//


#pragma once

#include "segmentation3d_base.h"

namespace misaxx_kidney_glomeruli {
    struct segmentation3d_klingberg : public segmentation3d_base {

        double m_max_glomerulus_radius = from_algorithm_json_or<double>("max-glomerulus-radius", 65);

        using segmentation3d_base::segmentation3d_base;

        void misa_work() override;
    };
}