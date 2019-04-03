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

#pragma once

#include "segmentation3d_base.h"

namespace misaxx_kidney_glomeruli {
    struct segmentation3d_klingberg : public segmentation3d_base {

        parameter<double> m_max_glomerulus_radius;

        using segmentation3d_base::segmentation3d_base;

        void work() override;

        void create_parameters(misaxx::misa_parameter_builder &t_parameters) override;

    };
}