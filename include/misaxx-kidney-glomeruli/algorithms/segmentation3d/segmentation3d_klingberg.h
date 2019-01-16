//
// Created by rgerst on 12.09.18.
//


#pragma once

#include <misaxx-kidney-glomeruli/algorithms/segmentation3d/segmentation3d_base.h>

namespace misaxx_kidney_glomeruli {
    struct segmentation3d_klingberg : public segmentation3d_base {

        parameter<double> m_max_glomerulus_radius;

        using segmentation3d_base::segmentation3d_base;

        void work() override;

        void create_parameters(misaxx::misa_parameter_builder &t_parameters) override;

    };
}