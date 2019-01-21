//
// Created by rgerst on 17.08.18.
//


#pragma once

#include <misaxx/core/misa_module.h>
#include "module_interface.h"

namespace misaxx_kidney_glomeruli {

    struct module : public misaxx::misa_module<module_interface> {
        using misaxx::misa_module<module_interface>::misa_module;

        parameter<std::string> m_segmentation2d_algorithm;
        parameter<std::string> m_segmentation3d_algorithm;
        parameter<std::string> m_quantification_algorithm;

        void create_blueprints(blueprint_list &t_blueprints, parameter_list &t_parameters) override;

        void build(const blueprint_builder &t_builder) override;

    };
}
