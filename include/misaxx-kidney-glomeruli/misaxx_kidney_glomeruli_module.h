//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/misa_module.h>
#include "../../src/misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_base.h"
#include "../../src/misaxx-kidney-glomeruli/algorithms/segmentation3d/segmentation3d_base.h"
#include "../../src/misaxx-kidney-glomeruli/algorithms/quantification/quantification_base.h"
#include "../../src/misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_klingberg.h"
#include "../../src/misaxx-kidney-glomeruli/algorithms/segmentation2d/segmentation2d_local_otsu.h"
#include "../../src/misaxx-kidney-glomeruli/algorithms/quantification/quantification_klingberg.h"
#include "../../src/misaxx-kidney-glomeruli/algorithms/quantification/quantification_constrained_klingberg.h"
#include "../../src/misaxx-kidney-glomeruli/algorithms/segmentation3d/segmentation3d_klingberg.h"
#include "kidney_glomeruli.h"

namespace misaxx_kidney_glomeruli {

    struct misaxx_kidney_glomeruli_module : public misaxx::misa_module<kidney_glomeruli> {
        using misaxx::misa_module<kidney_glomeruli>::misa_module;

        parameter<std::string> m_segmentation2d_algorithm;
        parameter<std::string> m_segmentation3d_algorithm;
        parameter<std::string> m_quantification_algorithm;

        void create_blueprints(blueprint_list &t_blueprints, parameter_list &t_parameters) override;

        void build(const blueprint_builder &t_builder) override;

    };
}
