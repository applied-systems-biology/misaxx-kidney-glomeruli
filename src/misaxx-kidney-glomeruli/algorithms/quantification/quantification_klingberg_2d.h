//
// Created by ruman on 23.03.19.
//

#pragma once

#include "quantification_base.h"

namespace misaxx_kidney_glomeruli {
    struct quantification_klingberg_2d : public quantification_base {

        parameter<double>  m_glomeruli_min_rad;
        parameter<double>  m_glomeruli_max_rad;

        using quantification_base::quantification_base;

        void work() override;

        void create_parameters(misaxx::misa_parameter_builder &t_parameters) override;
    };
}

