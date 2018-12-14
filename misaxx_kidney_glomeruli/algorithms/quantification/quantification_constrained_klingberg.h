//
// Created by rgerst on 17.09.18.
//


#pragma once

#include <coixx/objects/label_min_max_position.h>
#include "quantification_base.h"

namespace misaxx_kidney_glomeruli {
    struct quantification_constrained_klingberg : public quantification_base {

        double  m_glomeruli_min_rad = from_algorithm_json_or<double>("glomeruli-min-rad", 15);
        double  m_glomeruli_max_rad = from_algorithm_json_or<double>("glomeruli-max-rad", 65);

        using quantification_base::quantification_base;

        void misa_work() override;
    };
}